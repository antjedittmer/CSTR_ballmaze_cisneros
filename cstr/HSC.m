    function [htemp,stemp,Hess,ef,Cs,bamin,bamax] = HSC(x_aug_pred,hor,Q_,R_,P,x_k,xk,ref,lims,no_constraints)
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


    p1=[xk(1); x_aug_pred(1:nx:nx*hor)];
    p2=[xk(2); x_aug_pred(2:nx:nx*hor)];
    ref_N = ref(end-(nx-1):end);

    rho1 = p1(1);
    rho2 = p2(1);
            

     
    Ad = [-q/V-k0*exp(-E_R/rho2)         -(E_R/rho2^2)*k0*exp(-E_R/rho2)*rho1;
             -DH/(p*Cp)*k0*exp(-E_R/rho2)   -q/V-UA/(V*p*Cp)-(E_R/rho2^2)*DH/(p*Cp)*k0*exp(-E_R/rho2)*rho1];

    Beval = [zeros(nx/2,ni);
                            0;
             UA/(V*p*Cp*Ts)];

    Aeval = [zeros(nx/2) eye(nx/2);zeros(nx/2) Ad];
    
    % Euler
    A0 = Aeval*Ts+eye(nx);
    B0 = Ts*Beval;


    stemp = zeros(nx*hor,ni*hor);
    htemp = zeros(nx*hor,nx);        
    htemp(1:nx,1:nx) = A0;
    stemp(1:nx,1:ni) = B0;


    for mm=1:hor-1 
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
 
        stemp(mm*nx+1:(mm+1)*nx,1:mm*ni) = aeval*stemp((mm-1)*nx+1:mm*nx,1:mm*ni);     
        stemp(mm*nx+1:(mm+1)*nx,mm*ni+1:(mm+1)*ni) = beval;
        htemp(mm*nx+1:(mm+1)*nx,1:nx) = aeval*htemp((mm-1)*nx+1:mm*nx,:);
    end
        
    ctemp = kron(eye(hor),[1 0 0 0;0 1 0 0]);
    s_N = stemp(nx*(hor-1)+1:end,:);
    h_N = htemp(nx*(hor-1)+1:end,:);
    
    Hess = 2*(stemp'*Q_*stemp+R_+s_N'*P*s_N);
    ef = 2*((stemp'*Q_*htemp+s_N'*P*h_N)*x_k-(stemp'*Q_*ref+s_N'*P*ref_N));
    
    if no_constraints == 0
        Cs = ctemp*stemp;
        bamax = kron(ones(hor,1),[lims.Camax;lims.Tmax])-ctemp*htemp*x_k;
        bamin = kron(ones(hor,1),[lims.Camin;lims.Tmin])-ctemp*htemp*x_k;
    else
        Cs = [];
        bamin = [];
        bamax = [];
    end
    end
    