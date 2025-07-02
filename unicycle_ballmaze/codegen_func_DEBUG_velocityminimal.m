function [u,Xold,solv,xref,XX] = codegen_func_DEBUG_velocity(y,vel,iterations,u_lim,hor,x_aug_pred,Q_,R_,P,options,ref,Lc,u_old,C_delta,L,Cx,Cy)


x_k = [y;vel];
nx = 9;
solv = 0;
no = 0;
ni = 2;



%% MPC


if isempty(x_aug_pred)%zeros(120,1)
    x_aug_pred = kron(ones(hor,1),x_k);

else
    x_aug_pred(1:end-nx) = x_aug_pred(nx+1:end); 
end

% count = 0;
XX = [];
for ii = 1:iterations   
    %%% HSC
    [h,s,Hess,ef,Cs_,ba] = HSC_velocityminimal(x_aug_pred,hor,Q_,R_,P,x_k,ref,Lc,Cx,Cy);
    
    
    Cs = [C_delta;-C_delta;Cs_];
    uba = [u_lim-L*u_old;u_lim+L*u_old;ba];
%     [du,~,flag,~,~,auxout] = qpOASES(Hess,ef,Cs,[],[],[],uba,options);
    options = optimoptions('quadprog', 'Display', 'none');
    du = quadprog(Hess,ef,Cs,uba,[],[],[],[],[],options);

%     solv = solv + auxout.cpuTime;
    solv = solv + 0;

    if ~isempty(du)
        x_aug_pred =h*x_k+s*du;
    end
    XX = [XX x_aug_pred];
end

if ~isempty(du)
    u = (du(1:2)+u_old)';
else
    % if problem was infeasible, reuse last solution
    % Could be handled better by saving the whole u trajectory
    % and shifting it every infeasible timestep
    u = u_old';
end
Xold = x_aug_pred;
xref = ref;


end