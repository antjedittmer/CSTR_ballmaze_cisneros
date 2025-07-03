function [Lambda,S,Hess,ef,Cs,bamin,bamax] = HSC(X,N,Q_,R_,P,x_k,xk,ref,lims,no_constraints)
% HSC calculates matrices and boundaries for the QP
% otpimization for finding inputs for a q-LMPC with
% codegen_func_DEBUG_velocityminimal. This is for the CSTR example.
%
% Inputs
% - X: State prediction trajectory last sample
% - N: Samples prediction horizon
% - Q_,R_,P: Cost matrices
% - x_k: Current state
% - ref: Reference
% - lims: Constraints
% - no_constraints: Number of constraints (used as switch for
%   non-constrained optimization
%
% Output
% - the matrices Lambda and S for the state prediction trajectory,
% - the Hessian Hess and gradient ef
% - Boundaries bamin and bamax

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for presentation at DLR-FT Seminar: Antje Dittmer


nx = 4; ni=1;
Ts = 0.03;

%%%constants/nominal
q = 100;
V = 100;
p = 1000;
Cp = 0.239;
DH = -5e4;
E_R = 8750;
k0 = 7.2e10;
UA = 5e4;


p1=[xk(1); X(1:nx:nx*N)];
p2=[xk(2); X(2:nx:nx*N)];
ref_N = ref(end-(nx-1):end);

rho1 = p1(1);
rho2 = p2(1);


% CSTR model
Ad = [-q/V-k0*exp(-E_R/rho2)         -(E_R/rho2^2)*k0*exp(-E_R/rho2)*rho1;
    -DH/(p*Cp)*k0*exp(-E_R/rho2)   -q/V-UA/(V*p*Cp)-(E_R/rho2^2)*DH/(p*Cp)*k0*exp(-E_R/rho2)*rho1];

Beval = [zeros(nx/2,ni);
    0;
    UA/(V*p*Cp*Ts)];

Aeval = [zeros(nx/2) eye(nx/2);zeros(nx/2) Ad];

% Euler
A0 = Aeval*Ts+eye(nx);
B0 = Ts*Beval;

% Loop for building Lambda and S matrices
S = zeros(nx*N,ni*N);
Lambda = zeros(nx*N,nx);
Lambda(1:nx,1:nx) = A0;
S(1:nx,1:ni) = B0;


for mm=1:N-1
    rho1 = p1(mm+1);
    rho2 = p2(mm+1);

    Ad = [-q/V-k0*exp(-E_R/rho2)         -(E_R/rho2^2)*k0*exp(-E_R/rho2)*rho1;
        -DH/(p*Cp)*k0*exp(-E_R/rho2)   -q/V-UA/(V*p*Cp)-(E_R/rho2^2)*DH/(p*Cp)*k0*exp(-E_R/rho2)*rho1];

    Beval = [zeros(nx/2,ni);
        0;
        UA/(V*p*Cp*Ts)];

    Aeval = [zeros(nx/2) eye(nx/2);zeros(nx/2) Ad];

    % Euler
    aeval = Aeval*Ts+eye(nx);
    beval = Ts*Beval;

    S(mm*nx+1:(mm+1)*nx,1:mm*ni) = aeval*S((mm-1)*nx+1:mm*nx,1:mm*ni);
    S(mm*nx+1:(mm+1)*nx,mm*ni+1:(mm+1)*ni) = beval;
    Lambda(mm*nx+1:(mm+1)*nx,1:nx) = aeval*Lambda((mm-1)*nx+1:mm*nx,:);
end

ctemp = kron(eye(N),[1 0 0 0;0 1 0 0]);
s_N = S(nx*(N-1)+1:end,:);
h_N = Lambda(nx*(N-1)+1:end,:);

Hess = 2*(S'*Q_*S+R_+s_N'*P*s_N);
ef = 2*((S'*Q_*Lambda+s_N'*P*h_N)*x_k-(S'*Q_*ref+s_N'*P*ref_N));

if no_constraints == 0
    Cs = ctemp*S;
    bamax = kron(ones(N,1),[lims.Camax;lims.Tmax])-ctemp*Lambda*x_k;
    bamin = kron(ones(N,1),[lims.Camin;lims.Tmin])-ctemp*Lambda*x_k;
else
    Cs = [];
    bamin = [];
    bamax = [];
end
end
