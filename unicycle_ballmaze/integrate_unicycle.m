function xd = integrate_unicycle(t, x, u)
% integrate_unicycle models the unicycle dynamics.

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for presentation at DLR-FT Seminar: Antje Dittmer
    m = 1;
    I = 0.015;
    
    
    v = x(3);
    th = x(4);
    w = x(5);
    
    F = u(1);
    T = u(2);
    
    xd = zeros(5,1);
    xd(1) = v*cos(th);
    xd(2) = v*sin(th);
    xd(3) = 1/m*F;
    xd(4) = w;
    xd(5) = 0.075/I*T;

end