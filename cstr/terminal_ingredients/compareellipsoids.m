clear; clc; close all;

figDir = 'figDir';
if ~isfolder(figDir)
    mkdir(figDir)
end

ind1 = 1;
ind2 = 2;

load CSTR_FYZ;
W = @(rho) Wval{1}(rho);
Z = @(rho) Zval{1}(rho);

figure(1)

vol1 = [];
t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
z = [cos(t); sin(t)];

for i = 3:7 %3:length(paramgrid)
    for j = 4:8 %3:length(paramgrid2)
        Rchol = chol(W([paramgrid(i),paramgrid2(j)]));
        % Rchol = chol(W([0.5,350]));
        ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z; %#ok<*MINV> Small matrix
        plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'Color',[0.5 0.5 0.5],'LineWidth',1.5)
        hold on
        vol1 = [vol1 prod(svd(Z([paramgrid(i),paramgrid2(j)])))]; %#ok<AGROW> Small matrix
    end
end

volend= [];
W = @(rho) Wval{10}(rho);
Z = @(rho) Zval{10}(rho);
for i = 3:7%3:length(paramgrid)
    for j = 4:8%3:length(paramgrid2)
        Rchol = chol(W([paramgrid(i),paramgrid2(j)]));
        % Rchol = chol(W([0.5,350]));
        ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
        plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'Color','red','LineWidth',1.5)
        hold on
        volend = [volend prod(svd(Z([paramgrid(i),paramgrid2(j)])))];
    end
end

load CSTR_LMI;
vol2 = [];
for i = 3:7%:length(paramgrid)
    for j = 7:2:15%:length(paramgrid2)
        Rchol = chol(W([paramgrid(i),paramgrid2(j)]));
        % Rchol = chol(W([0.5,350]));
        ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
        plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'b','LineWidth',1.5)
        hold on
        vol2 = [vol2 prod(svd(Z([paramgrid(i),paramgrid2(j)])))];
    end
end


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

xlabel('C_A (mol/l)'); ylabel('T (K)'); grid on;

% black (const), grey (BMI 1), red (BM1 10), blue (LMI)
% tileSt = {'Constant ellipsoid', 'parameter dependent ellipsoids BMI approach 1st iteration', 'parameter dependent '...
% 'ellipsoids BMI approach 10th iteration','parameter dependent ellipsoid of the presented LMI approach'};

tileSt = {
    'Constant ellipsoid', ...
    'Parameter dependent ellipsoids BMI approach 1st iteration', ...
    'Parameter dependent ellipsoids BMI approach 10th iteration', ...
    'Parameter dependent ellipsoid of the presented LMI approach'
};

% Define colors (rows correspond to strings)
cl = [
    0     0     0;   % black
    0.5   0.5   0.5; % grey
    1     0     0;   % red
    0     0     1    % blue
];

% Helper function to convert RGB to LaTeX-style string
rgbStr = @(v) sprintf('%.2f,%.2f,%.2f', v);

% Build LaTeX string
titleCell = cell(numel(tileSt),1);
for i = 1:numel(tileSt)
    aStr = ['{\color[rgb]{', rgbStr(cl(i,:)), '}', tileSt{i}, '}'];
    titleCell{i} = aStr;
end
title(titleCell);

figDir = 'figDir';
if ~isfolder(figDir)
    mkdir(figDir)
end

strFig = 'Fig01_Ellipsoid';
print(fullfile(figDir,strFig), '-dpng');
%print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-depsc');

