ind1 = 1;
ind2 = 2;
clear clc

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH

load CSTR_FYZ
%W = @(rho) Wval{1}(rho);
%Z = @(rho) Zval{1}(rho);

W = @(rho) inv (Zval0 + basis {1} (rho) * Zval1 + basis {2} (rho) * Zval2);
Z = @(rho) (Zval0 + basis {1} (rho) * Zval1 + basis {2} (rho) * Zval2);
F = @(rho) (Fval0 + basis {1} (rho) * Fval1 + basis {2} (rho) * Fval2);
vol1= [];
for i = 3:7%3:length(paramgrid)
    for j = 4:8%3:length(paramgrid2)
Rchol = chol(W([paramgrid(i),paramgrid2(j)]));
% Rchol = chol(W([0.5,350]));
 t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
 z = [cos(t); sin(t)];
 ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
 plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'Color',[0.5 0.5 0.5],'LineWidth',1.5)
 hold on
 vol1 = [vol1 prod(svd(Z([paramgrid(i),paramgrid2(j)])))];
 
 
 x = 0.3:0.01:0.7;
f = F([paramgrid(i),paramgrid2(j)]);
plot(x,-f(1)/f(2)*(x-0.5)+20/f(2) + 350)
plot(x,-f(1)/f(2)*(x-0.5)-20/f(2) + 350)

    end
end

##volend= [];
##%W = @(rho) Wval{10}(rho);
##%Z = @(rho) Zval{10}(rho);
##W = @(rho) inv (Zval0 + basis {1} (rho) * Zval1 + basis {2} (rho) * Zval2);
##Z = @(rho) (Zval0 + basis {1} (rho) * Zval1 + basis {2} (rho) * Zval2);
##for i = 3:7%3:length(paramgrid)
##    for j = 4:8%3:length(paramgrid2)
##Rchol = chol(W([paramgrid(i),paramgrid2(j)]));
##% Rchol = chol(W([0.5,350]));
## t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
## z = [cos(t); sin(t)];
## ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
## plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'Color','red','LineWidth',1.5)
## hold on
## volend = [volend prod(svd(Z([paramgrid(i),paramgrid2(j)])))];
##    end
##end

##load CSTR_LMI
##vol2 = [];
##W = @(rho) a * inv (Y0 + basis {1} (rho) * Y1 + basis {2} (rho) * Y2);
##Z = @(rho) 1 / a * (Y0 + basis {1} (rho) * Y1 + basis {2} (rho) * Y2);
##for i = 3:7%:length(paramgrid)
##    for j = 7:2:15%:length(paramgrid2)
##Rchol = chol(W([paramgrid(i),paramgrid2(j)]));
##% Rchol = chol(W([0.5,350]));
## t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
## z = [cos(t); sin(t)];
## ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
## plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'b','LineWidth',1.5)
## hold on
## vol2 = [vol2 prod(svd(Z([paramgrid(i),paramgrid2(j)])))];
##    end
##end


load CSTR_LMI_constant
volc = [];
Rchol = chol(W(1));
% Rchol = chol(W([0.5,350]));
 t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
 z = [cos(t); sin(t)];
 ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
 plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'k--','LineWidth',1.5)
 hold on
 volc = [volc prod(svd(Z(1)))];

x = 0.3:0.01:0.7;
F = X0/Y0;
plot(x,-F(1)/F(2)*(x-0.5)+20/F(2) + 350)
plot(x,-F(1)/F(2)*(x-0.5)-20/F(2) + 350)
