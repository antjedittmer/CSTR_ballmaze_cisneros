function [u,Xold,solv,xref,XX] = codegen_func_DEBUG_velocityminimal(y,vel,iterations,u_lim,N,X,Q_,R_,P,optionsQ,ref,Lc,u_old,C_delta,L,Cx,Cy)
% codegen_func_DEBUG_velocityminimal calculates the optimal control input
% trajectory. 
% It calls HSC_velocityminimal to calculate 
% - the matrices Lambda and S for the state prediction trajectory, 
% - the Hessian Hess and gradient ef
% - Boundary matrix and vector Cs_ and ba
% These are then used to find the optimal control and state trajectory,
% either with qpOASES or quadprog.
% 
% Inputs
% - y: Output part of state vector x_k = [y;vel]
% - vel: Velocity part of state vector  x_k = [y;vel]
% - u_lim

x_k = [y;vel];
nx = 9;
solv = 0;
% no = 0;
% ni = 2;


%% MPC

if isempty(X)%zeros(120,1)
    X = kron(ones(N,1),x_k);

else
    X(1:end-nx) = X(nx+1:end);
end

% count = 0;
XX = [];
for ii = 1:iterations

    %%% HSC: Calculate 
    [Lambda,S,Hess,ef,Cs_,ba] = HSC_velocityminimal(X,N,Q_,R_,P,x_k,ref,Lc,Cx,Cy);
    Hess = 0.5*(Hess + Hess'); % Make sure Hessian is symmetric

    Cs = [C_delta; -C_delta; Cs_];
    uba = [u_lim-L*u_old;u_lim+L*u_old;ba];

    %     [du,~,flag,~,~,auxout] = qpOASES(Hess,ef,Cs,[],[],[],uba,optionsQ);
    % Options for quadprog
    options = optimoptions('quadprog', 'Display', 'none');

    du = quadprog(Hess,ef,Cs,uba,[],[],[],[],[],options);


    if ~isempty(du)
        X = Lambda*x_k + S*du;
    end
    XX = [XX X]; %#ok<*AGROW>
end

    %     solv = solv + auxout.cpuTime;
    solv = solv + 0;


if ~isempty(du)
    u = (du(1:2)+u_old)';
else
    % if problem was infeasible, reuse last solution
    % Could be handled better by saving the whole u trajectory
    % and shifting it every infeasible timestep
    u = u_old';
end
Xold = X;
xref = ref;

 end