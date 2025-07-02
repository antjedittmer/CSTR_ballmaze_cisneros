clc
clear


clc
clear


%% physical  params

Ts = 0.03;
horizon=25;
nx = 2;
ni = 1;

%% Pendulum on a cart

q = 100;
Tff = 350;
V = 100;
p = 1000;
Cp = 0.239;
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


%rho1 = t1, rho(2) = t1d

A_ = @(rho1,rho2) [-q/V-k0*exp(-E_R/rho2)         -(E_R/rho2^2)*k0*exp(-E_R/rho2)*rho1;
         -DH/(p*Cp)*k0*exp(-E_R/rho2)   -q/V-UA/(V*p*Cp)-(E_R/rho2^2)*DH/(p*Cp)*k0*exp(-E_R/rho2)*rho1];
B_ = @(rho1,rho2) [0;
         UA/(V*p*Cp*Ts)];

Blin = @(rho1,rho2) [0;
         UA/(V*p*Cp)];
     
%      A_ = @(rho) [0 [1 0];zeros(2,1) Ad(rho)];
% 
% B_ = @(rho) [0;0;
%          UA/(V*p*Cp*Ts)];

%      A_ = @(rho) [zeros(2) eye(2);zeros(2) Ad(rho)];
% 
% B_ = @(rho) [0;0;0;
%          UA/(V*p*Cp*Ts)];

    

A = @(rho1,rho2) eye(nx)+Ts*A_(rho1,rho2);

B2 = @(rho) Ts*B_(rho);

B = @(rho1,rho2) Ts*Blin(rho1,rho2);
%%
basis{1}=@(rho)rho(1);
basis{2}=@(rho)rho(1)/rho(2)^2;
% basis{1}=@(rho)rho(1);
% basis{2}=@(rho)rho(2);
% basis{3}=@(rho)rho(1)*sin(rho(2));
% basis{4}=@(rho)rho(1)*cos(rho(2));


no_rho=2;
paramgrid = 0.3:0.05:0.7;
paramgrid2 = 340:1:360;
paramgrid_d = -0.01:0.02:0.01;
paramgrid_d2 = -0.2:0.4:0.2;


Q = diag([1/0.5,1/350]);
% Q = diag([1,0.01,0.01]);
R = 1/300;
% Q = diag([1,.01,.01,.01]);
% R  = .01;

u_bar = 20;



%% Yalmip setup

yalmip('clear');
warning('off','YALMIP:strict')
SDPoptions = sdpsettings('savesolveroutput', 1, ...
                                 'savesolverinput' , 1, ...
                                 'verbose'         , 1, ...
                                 'solver'          ,'sdpt3');  

SDPoptions.sdpt3.maxit          = 150;

nx=size(A(1,1),1);
ni=size(B(1,1),2);


%% Find initial F that minimizez H2 norm


yalmip('clear');
Yk = sdpvar(nx,nx,'symmetric');
Yk1 = sdpvar(nx,nx,'symmetric');
Xk = sdpvar(ni,nx,'full');
X0 = sdpvar(ni,nx,'full');
X1 = sdpvar(ni,nx,'full');
X2 = sdpvar(ni,nx,'full');
X3 = sdpvar(ni,nx,'full');
Y0 = sdpvar(nx,nx,'symmetric');
Y1 = sdpvar(nx,nx,'symmetric');
Y2 = sdpvar(nx,nx,'symmetric');
Y3 = sdpvar(nx,nx,'symmetric');
t = sdpvar(1,1,'symmetric');
a = sdpvar(1,1,'symmetric');


Yk = Y0;
Xk = X0;
LMIconst = [];
for ii = 1:length(paramgrid)
    for jj = 1:length(paramgrid2)
%         for kk = 1:length(paramgrid_d)
%             for ll = 1:length(paramgrid_d)
                rho(1) = paramgrid(ii);     rho(2) = paramgrid2(jj);
%                 rhod(1) = paramgrid_d(kk);  rhod(2) = paramgrid_d2(ll);
                
%                 Yk  = Y0;% + basis{1}(rho)*Y1      + basis{2}(rho)*Y2;% + basis{3}(rho)*Y3;
%                 Yk1 = Y0 + basis{1}(rho)*Y1 + basis{2}(rho)*Y2;% + basis{3}(rho+rhod)*Y3;
%                 Xk    = X0;% + basis{1}(rho)*X1 + basis{2}(rho)*X2;
                
                ellip = [ a*u_bar^2     Xk;
                         Xk'            Yk];
              
                feasibility = [Yk                                          (A(paramgrid(ii),paramgrid2(jj))*Yk+B(paramgrid(ii),paramgrid2(jj))*Xk)' Yk          Xk';...
                   A(paramgrid(ii),paramgrid2(jj))*Yk+B(paramgrid(ii),paramgrid2(jj))*Xk   Yk                                                       zeros(nx)     zeros(nx,ni);...
                   Yk                                                    zeros(nx)                                                 inv(Q)             zeros(nx,ni);...
                   Xk                                                    zeros(ni,nx)                                              zeros(ni,nx)  inv(R)];    

               LMIconst = [LMIconst feasibility>0];
                             
%             end
%         end
    end
end
               %%%%% These are the lines to limit the size of the ellipsoid by the projection of rho onto x
                LMIconst = [LMIconst [1 0]*Yk*[1 0]'<=a*((paramgrid(end)-paramgrid(1))/2)^2];
                LMIconst = [LMIconst [0 1]*Yk*[0 1]'<=a*((paramgrid2(end)-paramgrid2(1))/2)^2];
               %%%%%%%%%

LMIconst = [LMIconst ellip>0 Yk>t a>0];
Optim = 0.1*a-100*t;
solution  = solvesdp(LMIconst,Optim, SDPoptions)


Y0=double(Y0);
% Y1=double(Y1);
% Y2=double(Y2);
% Yval3=double(Y3);
a = double(a);
Z0=Y0/a;
% Z1=Y1/a;
% Z2=Y2/a;
% Zval3=double(Z3);
X0=double(X0);
% X1=double(X1);
% X2=double(X2);


feas_YZ=isfeasible(LMIconst)

Y= @(rho)Y0;%(Y0+basis{1}(rho)*Y1+basis{2}(rho)*Y2);%+basis{3}(rho)*Yval3);
P= @(rho)inv(Y0);%+basis{1}(rho)*Y1+basis{2}(rho)*Y2);%+basis{3}(rho)*Yval3);
Z= @(rho)1/a*(Y0);%+basis{1}(rho)*Y1+basis{2}(rho)*Y2);%+basis{3}(rho)*Zval3);
W= @(rho)a*inv(Y0);%+basis{1}(rho)*Y1+basis{2}(rho)*Y2);%+basis{3}(rho)*Zval3);

%Plot projection of the ellipsoid
ind1 = 1;
ind2 = 2;

II = eye(2);
Zz = inv(W([0.5,350]));
Zproj = [II(:,ind1) II(:,ind2)]'*Zz*[II(:,ind1) II(:,ind2)];
Wproj = inv(Zproj);
Wch = chol(Wproj);
 ph = linspace(0, 2*pi, 100);
 z = [cos(ph);sin(ph)];
 ellipse = inv(Wch) * z;
 plot(ellipse(1,:), ellipse(2 ,:))

hold on
pause(0.01)
Wsvd=svd(W([0.5,350]));















