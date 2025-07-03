function [Atemp,Btemp] = ABC(rho1,rho2)
% ABC calculates the system and input matrix of the CSTR

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for presentation at DLR-FT Seminar: Antje Dittmer

%#codegen
%%%constants/nominal
q = 100;
V = 100;
p = 1000;
Cp = 0.239;
DH = -5e4;
E_R = 8750;
k0 = 7.2e10;
UA = 5e4;
nx = 4; ni=1;

Ts = 0.03;
Ad = [-q/V-k0*exp(-E_R/rho2)         -(E_R/rho2^2)*k0*exp(-E_R/rho2)*rho1;
    -DH/(p*Cp)*k0*exp(-E_R/rho2)   -q/V-UA/(V*p*Cp)-(E_R/rho2^2)*DH/(p*Cp)*k0*exp(-E_R/rho2)*rho1];

Beval = [zeros(nx/2,ni);
    0;
    UA/(V*p*Cp*Ts)];

Aeval = [zeros(nx/2) eye(nx/2);zeros(nx/2) Ad];


%             A2 = Aeval*Aeval;
%             A3 = Aeval*A2;
%             II = eye(nx);
%             Apre = (1.0/24.0)*Ts^4*A3+(1.0/6.0)*Ts^3*A2+0.5*Ts^2*Aeval+Ts*II;
%             A0 = Aeval*Apre+II;
%             B0 = Apre*Beval;
Atemp = Aeval*Ts+eye(nx);
Btemp = Ts*Beval;
end
