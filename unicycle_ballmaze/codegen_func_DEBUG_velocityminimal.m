function [u,Xold,solv,xref,XX] = codegen_func_DEBUG_velocityminimal(y,vel,iterations,u_lim,N,X,Q_,R_,P,ref,Lc,u_old,C_delta,L,Cx,Cy)
% codegen_func_DEBUG_velocityminimal calculates the optimal control input
% trajectory. It calls HSC_velocityminimal to calculate 
% - the matrices Lambda and S for the state prediction trajectory, 
% - the Hessian Hess and gradient ef
% - Boundary matrix and vector Cs_ and ba
% These are then used to find the optimal control and state trajectory,
% either with qpOASES or quadprog.
% 
% Inputs
% - y: Output part of state vector x_k = [y;vel]
% - vel: Velocity part of state vector  x_k = [y;vel]
% - iteration: number of iterations for prdiction calculation
% - u_lim: Limit on inputs
% - N: Samples prediction horizon
% - X: State prediction trajectory last sample
% - Q_,R_,P: Cost matrices
% - ref: Reference/steady state trajectories
% - Lc: Constant matrix used in constraint computations
% - uold: Input of last sample
% - C_delta: Part of constraint matrix:  Cs = [C_delta; -C_delta; Cs_]
% - L: Multiplier for constraint vector: [u_lim - L*u_old; u_lim + L*u_old; ba]
% - Cx, Cy: Obstacle centers
%
% Outputs
% - u: New control input
% - Xold: State prediction trajectory
% - solv: Solver time (currently unused)
% - xref: Reference trajectory
% - XX: State prediction trajectory for all iterations

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for presentation at DLR-FT Seminar: Antje Dittmer


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

    %%% HSC: Calculate Lambda, S, Hessian, gradient, boundaries
    [Lambda,S,Hess,ef,Cs_,ba] = HSC_velocityminimal(X,N,Q_,R_,P,x_k,ref,Lc,Cx,Cy);
    Hess = 0.5*(Hess + Hess'); % Make sure Hessian is symmetric

    Cs = [C_delta; -C_delta; Cs_];
    uba = [u_lim - L*u_old; u_lim + L*u_old; ba];

    %%% Calculate optimal trajectory
    % For now: Choose optimizer by simply commenting out the one not used
    % optionsQ = 0;
    % [du,~,flag,~,~,auxout] = qpOASES(Hess,ef,Cs,[],[],[],uba,optionsQ);
    
    options = optimoptions('quadprog', 'Display', 'none'); % Options for quadprog
    du = quadprog(Hess,ef,Cs,uba,[],[],[],[],[],options);

    if ~isempty(du) %If solution available: Update state trajectory
        X = Lambda*x_k + S*du;
    end
    XX = [XX X]; %#ok<*AGROW>
end

    %  solv = solv + auxout.cpuTime;
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