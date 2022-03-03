%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data-driven Distributed Operation of Electricity and Natural Gas Systems
% Christos ORDOUDIS, Viet Anh NGUYEN, Jalal KAZEMPOUR, Pierre PINSON, Daniel KUHN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main script
%

clear all; close all; clc;

% call startup to add the necessary path
startup;

tic

%%

% Data input
RTS_12node;

% DRO input

% Definition of value of \epsilon for individual chance constraints
% and the parameters for the ADMM algorithm 

Input_param.eps = 0.1; % \epsilon
Input_param.Max_iter = 10000; % MaxIter
Input_param.tolerance = 1e-2; % \eta
Input_param.rho = 0.001; % \rho
Input_param.restrict = 0.05; % restriction at isolated approach

% Number of individual runs (number of coupled datasets in the numerical
% study)

IR_max = 100;
IR_sim = 100;

% Number of out of sample data for each individual run (N') for testing
% dataset

OOS_max = 200;
OOS_sim = 100;

% Number of maximum sample size (N)

N_max = 1000;

% Number of sample data in training dataset (N)

N = 100;

% Total number of data 

Nscen = IR_max * (N_max + OOS_max);

% Generation of data 

rng(4,'twister')

% Number of wind farms

wf=[1:6];

% Loading the historical data for wind farms

wff=AV_AEMO2(:,wf);

% Cutting off very extreme values

cut_off_eps = 1e-2;
wff(wff<cut_off_eps) = cut_off_eps;
wff(wff>(1-cut_off_eps)) = 1 - cut_off_eps;

% Logit-normal transformation (Eq. (1) in ref. [31])

yy=log(wff./(1-wff));

% Calculation of mean and variance, note that we increase the mean to have
% higher wind penetration in our test-case

mu = mean(yy)+1.5;
sigma_m=cov(yy);
sigma_m=sigma_m./(std(yy)'*std(yy));

% Inverse of logit-normal transformation (Eq. (2) in ref. [31])

R = chol(sigma_m);
y = repmat(mu,Nscen,1) + randn(Nscen,size(WindDATA,1))*R;
Wind = (1+exp(-y)).^(-1);

% Checking correlation, mean and true mean of data

corrcoef(Wind);
mean(Wind);
true_mean_Wind = (1+exp(-mu)).^(-1);

% Reshaping the data structure

nWind = Wind';
nWind = reshape(nWind,size(WindDATA,1), N_max+OOS_max, IR_max);

% Initializing the matrices to gather final results

j = 1;

% For each coupled dataset, we pick N and N' samples
WPf_max = nWind(:,1:N_max,j)';
WPr_max = nWind(:,N_max+1:N_max+OOS_max,j)';
WPf = WPf_max(1:N,:);
WPr = WPr_max(1:OOS_sim,:);

    
% Build the corresponding data related to wind power production
all = [1:N];
system_info.Wscen = WPf(all,:)';
system_info.mu = mean(Wind)';%mean(wff)';%mean(system_info.Wscen,2); %system_info.xi = system_info.Wscen - repmat(system_info.mu, 1, size(system_info.Wscen,2));
system_info.cov = cov(Wind);%cov(wff);   

temp2 = wff' - repmat(system_info.mu, 1, size(wff,1));
temp = temp2';

% Build the corresponding data for RO
system_info.Wscen_RO = de2bi(0:2^size(system_info.Wmax,1)-1)';
system_info.Wexp_RO = mean(system_info.Wscen_RO,2); 
system_info.xi = system_info.Wscen_RO - repmat(system_info.mu, 1, size(system_info.Wscen_RO,2));
system_info.xi2 = [ones(1,size(system_info.xi,2));system_info.xi];
%si.exp_xi_xit = [1,si.mu';si.mu,si.cov+si.mu*si.mu'];
system_info.exp_xi_xit = [1,zeros(1,size(system_info.mu,1));zeros(size(system_info.mu,1),1),system_info.cov+zeros(size(system_info.mu,1),1)*zeros(1,size(system_info.mu,1))];
system_info.mu2 = [1;zeros(size(system_info.mu,1),1)];%[1;si.mu];

pl_wff = wff' - mean(wff)';
pl_wff = pl_wff';

pl_Wind = Wind' - mean(Wind)';
pl_Wind = pl_Wind';

Distr_wind;

[H, pValue, W] = swtest(pl_wff(1:5000,2),0.05);

display('DeCoupled_CC:');
DeCoup_EL_Gas_CC = DeCoup_CC(system_info, Input_param);
DeCoup_X_CC = DeCoup_EL_Gas_CC.X;
DeCoup_Y_CC = DeCoup_EL_Gas_CC.Y;
DeCoup_Q_CC = DeCoup_EL_Gas_CC.Q;
DeCoup_Obj_CC = DeCoup_EL_Gas_CC.Obj;

display('Coupled_CC:');
Coup_EL_Gas_CC = Coup_EL_Gas2_CC(system_info, Input_param);
Coup_X_CC = Coup_EL_Gas_CC.X;
Coup_Y_CC = Coup_EL_Gas_CC.Y;
Coup_Q_CC = Coup_EL_Gas_CC.Q;
Coup_Obj_CC = Coup_EL_Gas_CC.Obj;

display('ADMM_CC:');
ADMM_EL_Gas_CC = ADMM_CC(system_info, Input_param);
ADMM_X_CC = ADMM_EL_Gas_CC.X;
ADMM_Y_CC = ADMM_EL_Gas_CC.Y;
ADMM_Q_CC = ADMM_EL_Gas_CC.Q;
ADMM_Obj_CC = ADMM_EL_Gas_CC.Obj;
ADMM_pr_res_CC = ADMM_EL_Gas_CC.pr_res;
ADMM_dual_res_CC = ADMM_EL_Gas_CC.dual_res;
ADMM_iter_CC = ADMM_EL_Gas_CC.iter;
ADMM_Time_CC = ADMM_EL_Gas_CC.Time;


for i=1:ADMM_iter_CC
    plot_X_CC(i) = norm(ADMM_X_CC(:,:,i) - Coup_X_CC,2);
    plot_Y_CC(i) = norm(ADMM_Y_CC(:,:,i) - Coup_Y_CC,2);
    plot_Q_CC(i) = norm(ADMM_Q_CC(:,:,i) - Coup_Q_CC,2);
    plot_Obj_CC(i) = (ADMM_Obj_CC(i) - Coup_Obj_CC)/Coup_Obj_CC;
    plot_pr_res_CC(i) = norm(ADMM_pr_res_CC(:,:,i),2);
end

for i=1:ADMM_iter_CC
    plot_X_CC_fr(i) = norm(ADMM_X_CC(:,:,i) - Coup_X_CC,'fro');
    plot_Y_CC_fr(i) = norm(ADMM_Y_CC(:,:,i) - Coup_Y_CC,'fro');
    plot_Q_CC_fr(i) = norm(ADMM_Q_CC(:,:,i) - Coup_Q_CC,'fro');
    plot_pr_res_CC_fr(i) = norm(ADMM_pr_res_CC(:,:,i),'fro');
end

  
for j=1:ADMM_iter_CC-2
    plot_dual_res_CC(j) = norm(ADMM_dual_res_CC(:,:,j),2);
end

lw = 4;
mz = 15;
fs = 45;

figure(1)
plot(plot_X_CC,'-','LineWidth',lw,'Color',[0,0.45,0.74])
hold on
plot(plot_Y_CC,'-','LineWidth',lw,'Color',[0.96,0.3,0.33])
% hold on
% plot(plot_Q_CC,'-','LineWidth',lw,'Color',[0,0.6,0.2])
xlabel('Iteration number','Interpreter','latex','FontSize',fs);
ylabel('Convergence of system dispatch (MW - kcf)','Interpreter','latex','FontSize',fs);
set(gca,'ygrid','on')
set(gca, 'FontSize', fs)
legend({'$$\overline{\chi}$$','$$\overline{\upsilon}$$'},'Interpreter','latex','Location','north','Orientation','horizontal');
box on;

figure(2)
plot(plot_Obj_CC*100,'-','LineWidth',lw,'Color',[0,0.45,0.74])
xlabel('Iteration number','Interpreter','latex','FontSize',fs);
ylabel('Convergence of objective function (\%)','Interpreter','latex','FontSize',fs);
set(gca,'ygrid','on')
set(gca, 'FontSize', fs)
legend({'$$\overline{\gamma}$$'},'Interpreter','latex','Location','north','Orientation','horizontal');
box on;

figure(3)
plot(plot_pr_res_CC,'-','LineWidth',lw,'Color',[0,0.45,0.74])
xlabel('Iteration number','Interpreter','latex','FontSize',fs);
ylabel('Convergence of primal residual','Interpreter','latex','FontSize',fs);
set(gca,'ygrid','on')
set(gca, 'FontSize', fs)
legend({'$$\overline{\varrho}$$'},'Interpreter','latex','Location','north','Orientation','horizontal');
box on;

figure(4)
plot(plot_pr_res_CC,'-','LineWidth',lw,'Color',[0,0.45,0.74])
hold on 
plot(plot_dual_res_CC,'-','LineWidth',lw,'Color',[0.96,0.3,0.33])
xlabel('Iteration number','Interpreter','latex','FontSize',fs);
ylabel('Convergence of primal and dual residual','Interpreter','latex','FontSize',fs);
set(gca,'ygrid','on')
set(gca, 'FontSize', fs)
legend({'Primal residual', 'Dual residual'},'Location','north','Orientation','horizontal');
box on;


Time = toc





