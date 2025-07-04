clc; clear;

%% Code to generate Figure 2, closed-loop trajectorie, predictions and 
% ellipsoids for constant and parameter dependent case described in the 
% paper "Nonlinear Model Predictive Control for models in 
% quasi-LPV form"

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for presentation at DLR-FT Seminar: Antje Dittmer


%% Code to generate Figure 2, closed-loop trajectorie, predictions and ellipsoids for constant and parameter dependent case
figDir = 'figDir';
if ~isfolder(figDir)
    mkdir(figDir)
end

atk = 2; % timestep

figure(2)

%%  Load SimCT_LMIapproach_2it_nosaturatedU.mat
load SimCT_LMIapproach_2it_nosaturatedU.mat
plot(state_sim(atk:end,1),state_sim(atk:end,2),'--','Color',[0,0,0.8],'LineWidth',1.5)
hold on
plot(XX(1:4:end,atk),XX(2:4:end,atk),'Color',[0,0,0.8],'Marker','*','Markersize',6,'LineWidth',1.5)

load TerminalLMICT
ind1 = 1;
ind2 = 2;
Rchol = chol(W([.5,350]));
t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
z = [cos(t); sin(t)];
ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'Color',[0,0,0.8],'LineWidth',1.5)


%% Load SimLMI_2it_nosaturatedU
load SimLMI_2it_nosaturatedU

plot(state_sim(atk:end,1),state_sim(atk:end,2),'--','Color',[0.8,0,0],'LineWidth',1.5)
hold on
plot(XX(1:4:end,atk),XX(2:4:end,atk),'Color',[0.8,0,0],'Marker','o','Markersize',5,'LineWidth',1.5)

load TerminalLMI
ind1 = 1;
ind2 = 2;
Rchol = chol(W([XX(end-3,atk),XX(end-2,atk)]));
t = linspace(0, 2*pi, 100); % or any high number to make curve smooth
z = [cos(t); sin(t)];
ellipse = inv(Rchol([ind1,ind2],[ind1,ind2])) * z;
plot(ellipse(1,:)+0.5, ellipse(2,:)+350,'Color',[0.8,0,0],'LineWidth',1.5)


%% Add labels and title to plot
xlabel('C_A (mol/l)'); ylabel('T (K)'); grid on;

tileSt = {
    'Constant closed loop trajectories', ...
    'Constant prediction at k = 0', ...
    'Constant ellipsoid', ...
    'Parameter dependent LMIclosed loop trajectories', ...
    'Parameter dependent LMI prediction at k = 0', ...
    'Parameter dependent LMI approach ellipsoid'
};

fs  = 10.5;
legend(tileSt, 'FontSize',fs)


strFig = 'Fig02_ClLoop';
print(fullfile(figDir,strFig), '-dpng');