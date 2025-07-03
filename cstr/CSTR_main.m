clc;  clear; close all;

% Main function to calculate the CSTR qLMPC.

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH


X0 = [0.2 370];
Xref = [0.5 350];

Ts = 0.03;

% Add the Synthesis tool if necessary
if isempty(which('yalmip'))
    dirParent = fullfile(fileparts(pwd), 'SynthesisToolsP');
    addpath(genpath(dirParent));
end

%%%constants/nominal
q = 100;
Tff = 350;
V = 100;
rho = 1000;
crho = 0.239;
DH = -5e4;
E_R = 8750;
k0 = 7.2e10;
UA = 5e4;
Caf = 1;
Tcmin = 280;
Tcmax = 370;
Camin = 0;
Camax = 1;
Tmin = 280;
Tmax = 370;

lims.Tcmin = Tcmin;
lims.Tcmax = Tcmax;
lims.Camin = Camin-Xref(1);
lims.Camax = Camax-Xref(1);
lims.Tmin = Tmin-Xref(2);
lims.Tmax = Tmax-Xref(2);


%% MPC tuning

N = 10;
P = diag([100,100]);
Q = diag([1/0.5 1/350]);
R = 1/300;

Q = blkdiag(Q,zeros(2));
P = blkdiag(P,zeros(2));

Q_ = kron(eye(N),Q);
R_ = kron(eye(N),R);

% Set to true to use constant terminal ingredients (W, P)
use_constant_PW = false;

%% PARAMETERS SIMULATION


Xref = [Xref 0 0];
Xref = repmat(Xref,N,1);
ref = reshape(Xref',[4*N,1]);
ref = zeros(4*N,1);


%% SIMULATION LOOP
display('------------------------------------------------------------------')
display('               Simulation Loop'                                    )
display('------------------------------------------------------------------')

iter = 0; time = 0;
Tf = 100*Ts;
KKT_MPC = []; INFO_MPC = [];
controls_MPC = [];
state_sim = X0;

% visualize;
X = [];%zeros(4*N,1);
U = [];
ex_t = [];
sol_t=[];


options = [];
XX =[];
XXnl =[];
XXlpv = [];
w = [];

output.u = 0;
while time(end) < Tf
    tic
    if size(state_sim,1)>1
        vel = [(state_sim(end,1)-state_sim(end-1,1))/Ts;(state_sim(end,2)-state_sim(end-1,2))/Ts];
    else
        vel = zeros(2,1);
    end
    [output.u X solv Xlpv wout Xnl] = codegen_func_DEBUG(state_sim(end,1:2)',Xref(1,1:2)',vel,2,lims,N,X,Q_,R_,P,options,ref,output.u, use_constant_PW);
    XX = [XX X];
    XXnl = [XXnl Xnl];
    XXlpv{iter+1} = Xlpv;    
    w = [w wout];

    output.u(find(output.u>Tcmax)) = Tcmax;
    output.u(find(output.u<Tcmin)) = Tcmin;
    tt = toc;
    ex_t = [ex_t tt];

    % Simulate system
    sim_input.x = state_sim(end,:).';
    sim_input.u = output.u(1,:).';
    [~,xf] = ode45(@(t,x) integrate_cstr(t,x,sim_input.u),[0 Ts],sim_input.x);
    state_sim = [state_sim; xf(end,:)];
    
    iter = iter+1;
    nextTime = iter*Ts; 

    time = [time nextTime];
    U = [U; sim_input.u];
end

%% Save results
if use_constant_PW
    save('sim_CT.mat');
    disp('Saved simulation results in sim_CT.mat');
else
    save('sim_PD.mat');
    disp('Saved simulation results in sim_PD.mat'); 
end


%%
figure(1)
plot(state_sim(:,1))
figure(2)
plot(state_sim(:,2))
figure(3)
plot(U)
mean(ex_t)
