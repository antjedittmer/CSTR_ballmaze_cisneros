function [u,Xold,solv,Xlpv,w,Xnl] = codegen_func_DEBUG(y,yr,vel,iterations,lims,hor,x_aug_pred,Q_,R_,P,options,ref,u_old,use_constant_WP)

x_k = [y-yr;vel];
xk = [y;vel];
solv = 0;
no_term = 0;
ni = 1;
nx = 4;
C_delta = tril(kron(ones(hor),eye(ni)));
L = ones(hor,1);
%% P W

if use_constant_WP
    % constant: LMI approach (SimCT_LMIapproach_2it_nosaturatedU)
    W = 1e2*[1.864270491176637   0.035639680226742
       0.035639680226742   0.003130922854002];
    W = blkdiag(W,zeros(2));

    P =    1.0e+2 *[8.305160238354860   0.158771624894413
       0.158771624894413   0.013947984543811];
    P = blkdiag(P,zeros(2)); 
else
% LMI approach (SimLMI_2it_nosaturatedU.mat)
    Z0 =  [0.014145382964212  -0.450159879425466
      -0.450159879425466  24.225135696481786];

    Z1 =1e2*[-0.002868834775177   0.094968645936320
       0.094968645936320  -5.285023529129626];

    Z2 =      1.0e+07 *[0.003724071158792  -0.120290410004553
      -0.120290410004553   6.553988459860370];

    Y0 = [0.001795283989435  -0.057132763832769
      -0.057132763832769   3.074571990578813];
    Y1 = [-0.036410277143006   1.205309817197069
       1.205309817197069 -67.075724634936961];
    Y2 =    1.0e+06 *[0.004726464701457  -0.152668505129965
      -0.152668505129965   8.318099678669686];
end


%% MPC

if isempty(x_aug_pred)%zeros(120,1)
    x_aug_pred = kron(ones(hor,1),xk);
    no_term = 0;
else
    x_aug_pred(1:end-nx) = x_aug_pred(nx+1:end); 
end

Xlpv = [];
for ii = 1:iterations
    if ~use_constant_WP
        W = blkdiag(eye(2)/(Z0+x_aug_pred(end-3)*Z1+(x_aug_pred(end-3)/x_aug_pred(end-2)^2)*Z2),zeros(2));
        P = blkdiag(eye(2)/(Y0+x_aug_pred(end-3)*Y1+(1/x_aug_pred(end-2)^2)*Y2),zeros(2));
    end
    %%% HSC
    [h,s,Hess,ef,Cs_,bamin,bamax] = HSC(x_aug_pred,hor,Q_,R_,P,x_k,xk,ref,lims,no_term);

    s_N = s(nx*(hor-1)+1:end,:);
    h_N = h(nx*(hor-1)+1:end,:);    
    yalmip('clear')
    du = sdpvar(ni*hor,1);

    if no_term == 1; %dont use terminal constraint on first iteration
        diagnostic=optimize([C_delta*du+L*u_old;-C_delta*du-L*u_old]<=[L*lims.Tcmax;L*lims.Tcmin],...
        (h*x_k+s*du)'*Q_*(h*x_k+s*du)+du'*R_*du + (h_N*x_k+s_N*du)'*P*(h_N*x_k+s_N*du),sdpsettings('solver','sdpt3','verbose',0))
    else
        diagnostic = optimize([[C_delta*du+L*u_old;-C_delta*du-L*u_old]<=[L*lims.Tcmax;L*lims.Tcmin] (h_N*x_k+s_N*du)'*W*(h_N*x_k+s_N*du)<=1],...
        (h*x_k+s*du)'*Q_*(h*x_k+s*du)+du'*R_*du + (h_N*x_k+s_N*du)'*P*(h_N*x_k+s_N*du),sdpsettings('solver','sdpt3','verbose',0))  
    end
    du = double(du);

    du(find(du(~isfinite(du)))) = 0;

    x_aug_pred = h*(xk)+s*du;
    
    % lpv prediction
    Xnl = zeros(hor+1,4);
    Xnl(1,:) = xk';
    for jj = 2:hor+1
        [aeval,beval] = ABC(Xnl(jj-1,1),Xnl(jj-1,2));
        Xnl(jj,:) =  (aeval*Xnl(jj-1,:)'+beval*du(jj-1))';%+[yr;0;0]';
    end
    Xnl = reshape(Xnl(2:end,:)',[hor*4 1]);
    Xlpv = [Xlpv x_aug_pred];
    
    %lpv prediction
    no_term = 0;

end
u = du(1)+u_old;
Xold = x_aug_pred;
w = (h_N*x_k+s_N*du)'*W*(h_N*x_k+s_N*du);

end