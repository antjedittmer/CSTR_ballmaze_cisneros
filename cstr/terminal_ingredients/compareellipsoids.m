clear; clc; close all;

%% Code to generate Figure 1, ellipsoids generated with parameter dependent BMI, LIMI and constant LMI
% CSTR_FYZ.mat, CSTR_LMI.mat, and CSTR_LMI_constant.mat ara included in the
% git directory, but can me reproduces with the respective m-files

figDir = 'figDir';
if ~isfolder(figDir)
    mkdir(figDir)
end

ind1 = 1;
ind2 = 2;

vol1 = [];
t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
z = [cos(t); sin(t)];

figure(1)

%% Load CSTR_FYZ (can be generated with F_YZ.m): BMI approach, grey line and red line 
load CSTR_FYZ;
W = @(rho) Wval{1}(rho);
Z = @(rho) Zval{1}(rho);

for i = 3:7 % 3:length(paramgrid)
    for j = 4:8 % 3:length(paramgrid2)
        Rchol = chol(W([paramgrid(i),paramgrid2(j)]));
        % Rchol = chol(W([0.5,350]));
        ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z; %#ok<*MINV> Small matrix, non efficient inv is ok
        plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'Color',[0.5 0.5 0.5],'LineWidth',1.5) % grey ellipse
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
        plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'Color','red','LineWidth',1.5) % red ellipse
        hold on
        volend = [volend prod(svd(Z([paramgrid(i),paramgrid2(j)])))];
    end
end


%% Load CSTR_LMI (can be generated with LMIapproach.m): parameter dependent LMI approach, blue line
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

%% Load CSTR_LMI_constant (can be generated with LMIapproach_constant.m): LMI approach, blue line
load CSTR_LMI_constant
volc = [];
Rchol = chol(W(1));
% Rchol = chol(W([0.5,350]));
ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'k--','LineWidth',1.5)
hold on
volc = [volc prod(svd(Z(1)))];

%% Add labels and title to block

xlabel('C_A (mol/l)'); ylabel('T (K)'); grid on;

% Define the title cell
titleSt = {
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

% Helper function to convert RGB to string
rgbStr = @(v) sprintf('%.2f,%.2f,%.2f', v);

% Build LaTeX string
titleCell = cell(numel(titleSt),1);
for i = 1:numel(titleSt)
    aStr = ['{\color[rgb]{', rgbStr(cl(i,:)), '}', titleSt{i}, '}'];
    titleCell{i} = aStr;
end
title(titleCell);

figDir = 'figDir';
if ~isfolder(figDir)
    mkdir(figDir)
end

strFig = 'Fig01_Ellipsoid';
print(fullfile(figDir,strFig), '-dpng');
% print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-depsc');

