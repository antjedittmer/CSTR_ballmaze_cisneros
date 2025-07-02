function xd = integrate_cstr(t, x, u) %#ok<INUSD>

q = 100;
Tf = 350;
V = 100;
rho = 1000;
crho = 0.239;
DH = -5e4;
E_R = 8750;
k0 = 7.2e10;
UA = 5e4;
Caf = 1;

Ca = x(1);
T = x(2);

Tc = u(1);

xd = zeros(2,1);
xd(1) = q/V*(Caf-Ca)-k0*exp(-E_R/T)*Ca;
xd(2) = q/V*(Tf-T)+UA/(rho*V*crho)*(Tc-T)-DH*k0/(rho*crho)*exp(-E_R/T)*Ca;

end