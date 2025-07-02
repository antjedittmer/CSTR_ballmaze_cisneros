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
paramgrid2 = 340:2:360;
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
Y0 = sdpvar(nx,nx,'symmetric');
Y1 = sdpvar(nx,nx,'symmetric');
Y2 = sdpvar(nx,nx,'symmetric');
Xk = sdpvar(ni,nx,'full');
X0 = sdpvar(ni,nx,'full');
X1 = sdpvar(ni,nx,'full');
X2 = sdpvar(ni,nx,'full');
X3 = sdpvar(ni,nx,'full');
t = sdpvar(1,1,'symmetric');


LMIconst = [];
count = 1;
for ii = 1:length(paramgrid)
    for jj =1:length(paramgrid2)
        for kk = 1:length(paramgrid_d)
            for ll = 1:length(paramgrid_d)
                rh1   = paramgrid(ii); rh2 = paramgrid2(jj); rh = [rh1 rh2];
                rhd1  = paramgrid_d(kk); rhd2 = (paramgrid_d(ll)); rhd = [rhd1 rhd2];
                Yk    = Y0 + basis{1}(rh)*Y1 + basis{2}(rh)*Y2;
                Yk1   = Y0 + basis{1}(rh+rhd)*Y1 + basis{2}(rh+rhd)*Y2;
                Xk    = X0 + basis{1}(rh)*X1 + basis{2}(rh)*X2;% + basis{3}(rh)*X3;

                feasibility = [Yk                                          (A(paramgrid(ii),paramgrid2(jj))*Yk+B(paramgrid(ii),paramgrid2(jj))*Xk)' Yk          Xk';...
                    A(paramgrid(ii),paramgrid2(jj))*Yk+B(paramgrid(ii),paramgrid2(jj))*Xk   Yk1                                                       zeros(nx)     zeros(nx,ni);...
                    Yk                                                    zeros(nx)                                                 inv(Q)             zeros(nx,ni);...
                    Xk                                                    zeros(ni,nx)                                              zeros(ni,nx)  inv(R)];

                LMIconst = [LMIconst feasibility>0 Yk>=t*eye(nx)];
            end
        end
    end
end

Optim = -t;
solution  = solvesdp(LMIconst,[], SDPoptions)

Y0=double(Y0);
Y1=double(Y1);
Y2=double(Y2);
X0=double(X0);
X1=double(X1);
X2=double(X2);
% X3=double(X3)


Fval{1} = @(rho)(X0+basis{1}(rho)*X1+basis{2}(rho)*X2)/(Y0+basis{1}(rho)*Y1+basis{2}(rho)*Y2)
Pval{1} = @(rho) eye(nx)/(Y0+basis{1}(rho)*Y1+basis{2}(rho)*Y2)
tpd = double(t)
tval{1} = 0;
feas_F(1) = isfeasible(LMIconst)

%number of iterations
iterations = 20;
%%
for it = 1:iterations
    % fix F, find new Y and Z, optimize logdet(Z)

    F = @(rho)Fval{it}(rho);
    yalmip('clear');
    Z0 = sdpvar(nx,nx,'symmetric');
    Z1 = sdpvar(nx,nx,'symmetric');
    Z2 = sdpvar(nx,nx,'symmetric');
    Z3 = sdpvar(nx,nx,'symmetric');
    Y0 = sdpvar(nx,nx,'symmetric');
    Y1 = sdpvar(nx,nx,'symmetric');
    Y2 = sdpvar(nx,nx,'symmetric');
    Y3 = sdpvar(nx,nx,'symmetric');
    t = sdpvar(1,1,'symmetric');



    LMIconst = [];
    for ii = 1:length(paramgrid)
        for jj = 1:length(paramgrid2)
            for kk = 1:length(paramgrid_d)
                for ll = 1:length(paramgrid_d)
                    rho(1) = paramgrid(ii);     rho(2) = paramgrid2(jj);
                    rhod(1) = paramgrid_d(kk);  rhod(2) = paramgrid_d2(ll);

                    Yk  = Y0 + basis{1}(rho)*Y1      + basis{2}(rho)*Y2;% + basis{3}(rho)*Y3;
                    Yk1 = Y0 + basis{1}(rho+rhod)*Y1 + basis{2}(rho+rhod)*Y2;% + basis{3}(rho+rhod)*Y3;
                    Zk  = Z0 + basis{1}(rho)*Z1      + basis{2}(rho)*Z2;%  + basis{3}(rho)*Z3;
                    Zk1 = Z0 + basis{1}(rho+rhod)*Z1 + basis{2}(rho+rhod)*Z2;% + basis{3}(rho+rhod)*Z3;

                    ellip = [ Zk                                     Zk*(A(paramgrid(ii),paramgrid2(jj))+B(paramgrid(ii),paramgrid2(jj))*F(rho))';
                        (A(paramgrid(ii),paramgrid2(jj))+B(paramgrid(ii),paramgrid2(jj))*F(rho))*Zk    Zk1];

                    feasibility = [Yk                                                          Yk*(A(paramgrid(ii),paramgrid2(jj))+B(paramgrid(ii),paramgrid2(jj))*F(rho))'  Yk             Yk*F(rho)';...
                        (A(paramgrid(ii),paramgrid2(jj))+B(paramgrid(ii),paramgrid2(jj))*F(rho))*Yk   Yk1                                                            zeros(nx)     zeros(nx,ni);...
                        Yk                                                          zeros(nx)                                                    inv(Q)        zeros(nx,ni);...
                        F(rho)*Yk                                                    zeros(ni,nx)                                                 zeros(ni,nx)  inv(R)];

                    LMIconst = [LMIconst feasibility>0 ellip>=0 geomean(Zk)>t Zk>0];

                    %%%%% These are the lines to limit the size of the ellipsoid by the projection of rho onto x
                    LMIconst = [LMIconst [1 0]*Zk*[1 0]'<=((paramgrid(end)-paramgrid(1))/2)^2];
                    LMIconst = [LMIconst [0 1]*Zk*[0 1]'<=((paramgrid2(end)-paramgrid2(1))/2)^2];
                    %%%%%%%%%
                    LMIconst = [LMIconst ellip>=0 geomean(Zk)>t];

                    for mm = 1:size(F,1)
                        Fij = F(rho);
                        LMIconst = [LMIconst Fij(mm,:)*Zk*Fij(mm,:)'<=u_bar(mm)^2];
                    end
                end
            end
        end
    end


    LMIconst = [LMIconst t>0];
    Optim = -t;
    solution  = solvesdp(LMIconst,Optim, SDPoptions);


    Yval0=double(Y0);
    Yval1=double(Y1);
    Yval2=double(Y2);
    % Yval3=double(Y3);
    Zval0=double(Z0);
    Zval1=double(Z1);
    Zval2=double(Z2);
    % Zval3=double(Z3);


    feas_YZ(it)=isfeasible(LMIconst)

    Yval{it}= @(rho)(Yval0+basis{1}(rho)*Yval1+basis{2}(rho)*Yval2);%+basis{3}(rho)*Yval3);
    Pval{it}= @(rho)inv(Yval0+basis{1}(rho)*Yval1+basis{2}(rho)*Yval2);%+basis{3}(rho)*Yval3);
    Zval{it}= @(rho)(Zval0+basis{1}(rho)*Zval1+basis{2}(rho)*Zval2);%+basis{3}(rho)*Zval3);
    Wval{it}= @(rho)inv(Zval0+basis{1}(rho)*Zval1+basis{2}(rho)*Zval2);%+basis{3}(rho)*Zval3);

    %Plot projection of the ellipsoid
    ind1 = 1;
    ind2 = 2;

    II = eye(2);
    Z = inv(Wval{it}([0.5,350]));
    Zproj = [II(:,ind1) II(:,ind2)]'*Z*[II(:,ind1) II(:,ind2)];
    Wproj = inv(Zproj);
    Wch = chol(Wproj);
    ph = linspace(0, 2*pi, 100);
    z = [cos(ph);sin(ph)];
    ellipse = inv(Wch) * z;
    plot(ellipse(1,:), ellipse(2 ,:))

    hold on
    pause(0.01)
    Wsvd{it}=svd(Wval{it}([0.5,350]));


    %% fix Y,Z, find F optimize t-s

    yalmip('clear');
    F0 = sdpvar(ni,nx,'full');
    F1 = sdpvar(ni,nx,'full');
    F2 = sdpvar(ni,nx,'full');
    F3 = sdpvar(ni,nx,'full');
    t = sdpvar(1,1,'symmetric');
    s = sdpvar(1,1,'symmetric');


    LMIconst = [];
    for ii = 1:length(paramgrid)
        for jj = 1:length(paramgrid2)
            for kk = 1:length(paramgrid_d)
                for ll = 1:length(paramgrid_d)
                    rho(1) = paramgrid(ii); rho(2) = paramgrid2(jj);
                    rhod(1) = paramgrid_d(kk); rhod(2) = paramgrid_d2(ll);

                    Zk  = Zval{it}(rho);
                    Zk1 = Zval{it}(rho+rhod);
                    Yk  = Yval{it}(rho);
                    Yk1 = Yval{it}(rho+rhod);
                    F = (F0+basis{1}(rho)*F1+basis{2}(rho)*F2);%+basis{3}(rho)*F3);%/Yval{it}(rho);%+basis{3}(rho)*F3;


                    ellip = [ Zk                                     Zk*(A(paramgrid(ii),paramgrid2(jj))+B(paramgrid(ii),paramgrid2(jj))*F)';
                        (A(paramgrid(ii),paramgrid2(jj))+B(paramgrid(ii),paramgrid2(jj))*F)*Zk    Zk1];

                    feasibility = [Yk                                                          Yk*(A(paramgrid(ii),paramgrid2(jj))+B(paramgrid(ii),paramgrid2(jj))*F)'  Yk             Yk*F';...
                        (A(paramgrid(ii),paramgrid2(jj))+B(paramgrid(ii),paramgrid2(jj))*F)*Yk   Yk1                                                            zeros(nx)     zeros(nx,ni);...
                        Yk                                                          zeros(nx)                                                    inv(Q)        zeros(nx,ni);...
                        F*Yk                                                    zeros(ni,nx)                                                 zeros(ni,nx)  inv(R)];

                    LMIconst = [LMIconst feasibility>0 ellip>=0];
                    for mm = 1:size(F,1)
                        saturation = [u_bar(mm)^2+t  F(mm,:)*Zk;
                            Zk*F(mm,:)'          Zk];
                        LMIconst = [LMIconst saturation>=0];
                    end
                end
            end
        end
    end


    LMIconst = [LMIconst t<0];

    Optim = t;
    solution  = solvesdp(LMIconst,Optim, SDPoptions);

    Fval0 = double(F0);
    Fval1 = double(F1);
    Fval2 = double(F2);
    Fval3=double(F3);

    feas_F(it+1) = isfeasible(LMIconst)
    tval{it+1} = double(t)

    Fval{it+1}  = @(rho)(Fval0+basis{1}(rho)*Fval1+basis{2}(rho)*Fval2);%+basis{3}(rho)*Fval3);%/(Yval0+basis{1}(rho)*Yval1+basis{2}(rho)*Yval2);%+basis{3}(rho)*Fval3);

end


save('CSTR_FYZ'); 












