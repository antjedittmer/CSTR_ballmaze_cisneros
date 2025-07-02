clc;
clear all;
close all;

% Input limits
u_lim = 1;

% Obstacle centers
Cx = [1/3-0.5, 1/3-0.5, 1/3-0.5,   1-0.5, 1-0.5,   1-0.5, 5/3-0.5, 5/3-0.5, 5/3-0.5];
Cy = [1/3-1,   1-1, 5/3-1, 1/3-1, 1-1, 5/3-1, 1/3-1,   1-1, 5/3-1];
% Obstacle radii
radius = 0.3;

% Initial condition
X0 = [0.16 0.3 0 pi/2 0];

% Sampling time
Ts = 0.025;

% Physical parameters
m = 1;
I = 0.015;
l = 0.1;

%% MPC tuning

N = 30;
P = diag([1000 1000 .1 1 0 0 0 0 0]);
Q = diag([100 100 .01 .1 0.01 0.01 1 0.01 0.01]); 
R  = diag([.01, .1]);

Q_ = kron(eye(N),Q);
R_ = kron(eye(N),R);

%% SIMULATION LOOP
display('------------------------------------------------------------------')
display('               Simulation Loop'                                    )
display('------------------------------------------------------------------')

iter = 0; time = 0;
Tf = 75*10.5*Ts;
state_sim = X0;

% visualize;
X = [];
U = [];
ex_t = [];
sol_t=[];

ulim = u_lim*ones(2*N,1);
Xref = [0 0 0 0 0 0 0 0 0];
Xref = repmat(Xref,N,1);
ref = reshape(Xref',[9*N,1]);
%%%
reftime = 0:Ts:22*Tf;
refx = sin(1/3*reftime)';
refy = sin(2/3*reftime)';
refv = zeros(size(refx));
REF = zeros(9*(22*Tf/Ts+1),1);
REF(1:9:end) = refx;
REF(2:9:end) = refy;
%%%
% precompute constant matrix used in constraint computations
Lc = [];
for iiii = 1:length(Cx)
    Lc = [Lc;-radius^2 + Cx(iiii)^2 + Cy(iiii)^2]; %circle
end
Lc = kron(ones(N,1),Lc);
XX =[];
XXlpv = [];
XXref = [];

output.u = [0 0];
C_delta = tril(kron(ones(N),eye(2)));
L = kron(ones(N,1),eye(2));
while time(end) < Tf
    ref = REF(9*iter+1:9*iter+N*9);
    tic
    if size(state_sim,1)>1
        vel = [(state_sim(end,1)-state_sim(end-1,1))/Ts;(state_sim(end,2)-state_sim(end-1,2))/Ts;(state_sim(end,3)-state_sim(end-1,3))/Ts;state_sim(end,5);(state_sim(end,5)-state_sim(end-1,5))/Ts];
    else
        vel = zeros(5,1);
    end
    [output.u X solv xxref Xit] = codegen_func_DEBUG_velocityminimal(state_sim(end,1:4)',vel,2,ulim,N,X,Q_,R_,P,0,ref,Lc,output.u',C_delta,L,Cx,Cy);
    tt = toc
    XX = [XX X];
%     XXit{iter+1} = Xit;
    XXref = [XXref xxref];
    sol_t = [sol_t solv];
    
    
    if output.u(1) > u_lim
        output.u(1) = u_lim;
    end
    if output.u(1) < -u_lim;
        output.u(1) = -u_lim;
    end
    if output.u(2) > u_lim
        output.u(2) = u_lim;
    end
    if output.u(2) < -u_lim;
        output.u(2) = -u_lim;
    end

    ex_t = [ex_t tt];
    
    % Simulate system
    sim_input.x = state_sim(end,:).';
    sim_input.u = output.u(1,:).';
    [~,xf] = ode45(@(t,x) integrate_unicycle(t,x,sim_input.u),[0 Ts],sim_input.x);
    state_sim = [state_sim; xf(end,:)];
    
    iter = iter+1;
    nextTime = iter*Ts; 

    time = [time nextTime];
    
    U = [U; sim_input.u];
end
%% time plots

figure(1)
plot(state_sim(:,1))
hold on
plot(state_sim(:,2))
figure(2)
plot(U(1:2:end))
hold on
plot(U(2:2:end))
ylim([-1.5 1.5])
mean(ex_t)

%% x-y plots
figure()
animate_plot = false;

plotmaze_ref