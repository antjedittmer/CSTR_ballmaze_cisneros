function [Lambda,S,Hess,ef,Cs,ba] = HSC_velocityminimal(X,N,Q_,R_,P,x_k,ref,Lc,Cx,Cy)
% HSC_velocityminimal calculates matrices and boundaries for the QP
% otpimization for finding inputs for a q-LMPC with
% codegen_func_DEBUG_velocityminimal. This is for the unicycle example.
%
% Inputs
% - X: State prediction trajectory last sample
% - N: Samples prediction horizon
% - Q_,R_,P: Cost matrices
% - x_k: Current state
% - ref: Reference
% - uold: Input of last sample
% - C_delta: Part of constraint matrix:  Cs = [C_delta; -C_delta; Cs_]
% - L: Multiplier for constraint vector: [u_lim - L*u_old; u_lim + L*u_old; ba]
% - Cx, Cy: Obstacle centers
%
% Output
% - the matrices Lambda and S for the state prediction trajectory, 
% - the Hessian Hess and gradient ef
% - Boundary matrix and vector Cs_ and ba

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for presentation at DLR-FT Seminar: Antje Dittmer

% number of states and input
nx = 9; ni = 2;

% Values for Cx, Cy
% Cx = [1/3-0.5, 1/3-0.5, 1/3-0.5,   1-0.5, 1-0.5,   1-0.5, 5/3-0.5, 5/3-0.5, 5/3-0.5];
% Cy = [1/3-1,   1-1, 5/3-1, 1/3-1, 1-1, 5/3-1, 1/3-1,   1-1, 5/3-1];
Ts = 0.025;
ncirc = length(Cx);

% parameter trajectories
p1 = [x_k(1); X(1:nx:nx*N)];
p2 = [x_k(2); X(2:nx:nx*N)];
t1 = [x_k(4); X(4:nx:nx*N)];
t2 = [x_k(3); X(3:nx:nx*N)];
ref_N = ref(end-(nx-1):end);

rho1 = p1(2);
rho2 = p2(2);
rhot = t1(1);
rhov = t2(1);

% Unicycle model
Ad = [0 0 cos(rhot)  -rhov*sin(rhot) 0;
    0 0 sin(rhot)  rhov*cos(rhot)  0;
    0 0 0          0               0;
    0 0 0          0               1;
    0 0 0          0               0];

O = zeros(4);
I_ = eye(5);
I_ = I_(1:4,1:5);

Aeval = [O I_; zeros(5,4) Ad];
Beval = [zeros(4,ni);
    0 , 0;
    0 , 0;
    1/Ts , 0;
    0 , 0;
    0 , 0.075/(Ts*0.015)];


% 4th order discretization
A2 = Aeval*Aeval;
A3 = Aeval*A2;
II = eye(nx);
Apre = (1.0/24.0)*Ts^4*A3+(1.0/6.0)*Ts^3*A2+0.5*Ts^2*Aeval+Ts*II;
A0 = Aeval*Apre+II;
B0 = Apre*Beval;

C0 =[];
baeval = [];
for circ = 1:ncirc
    C0 = [C0;-2*rho1+2*Cx(circ), -2*rho2+2*Cy(circ), 0, 0,0,0,0,0,0]; %#ok<*AGROW> %circle *circ
    baeval = [baeval; -rho1^2+2*Cx(circ)*rho1-rho2^2+2*Cy(circ)*rho2];
end

% Loop for building Lambda and S matrices
ctemp = zeros(N*ncirc,nx*N);
S = zeros(nx*N,ni*N);
Lambda = zeros(nx*N,nx);
Lambda(1:nx,1:nx) = A0;
S(1:nx,1:ni) = B0;
ctemp(1:ncirc,1:nx) = C0;

for mm = 1:N-1
    rho1 = p1(mm+2);
    rho2 = p2(mm+2);
    rhot = t1(mm+1);
    rhov = t2(mm+1);

    Ad = [0 0 cos(rhot)  -rhov*sin(rhot) 0;
        0 0 sin(rhot)  rhov*cos(rhot)  0;
        0 0 0          0               0;
        0 0 0          0               1;
        0 0 0          0               0];

    Aeval = [O I_; zeros(5,4) Ad];

    A2 = Aeval*Aeval;
    A3 = Aeval*A2;
    II = eye(nx);
    Apre = (1.0/24.0)*Ts^4*A3+(1.0/6.0)*Ts^3*A2+0.5*Ts^2*Aeval+Ts*II;
    aeval = Aeval*Apre+II;
    beval = Apre*Beval;

    ceval = [];
    for circ = 1:ncirc
        ceval = [ceval;-2*rho1+2*Cx(circ), -2*rho2+2*Cy(circ), 0, 0,0,0,0,0,0];
        baeval = [baeval; -rho1^2+2*Cx(circ)*rho1-rho2^2+2*Cy(circ)*rho2];
    end

    S(mm*nx+1:(mm+1)*nx,1:mm*ni) = aeval*S((mm-1)*nx+1:mm*nx,1:mm*ni);
    S(mm*nx+1:(mm+1)*nx,mm*ni+1:(mm+1)*ni) = beval;
    Lambda(mm*nx+1:(mm+1)*nx,1:nx) = aeval*Lambda((mm-1)*nx+1:mm*nx,:);
    ctemp(mm*ncirc+1:mm*ncirc+ncirc,mm*nx+1:(mm+1)*nx) = ceval;
end

s_N = S(nx*(N-1)+1:end,:);
h_N = Lambda(nx*(N-1)+1:end,:);

Hess = 2*(S'*Q_*S+R_ + s_N'*P*s_N);
ef = 2*((S'*Q_*Lambda+s_N'*P*h_N)*x_k-(S'*Q_*ref+s_N'*P*ref_N));

Cs = ctemp*S;
ba = Lc + ctemp*X-ctemp*Lambda*x_k-baeval; %the additional ctemp*x_aug_pred correspondons to -grad(f)*rho from grad(f)*(x-rho)
end
