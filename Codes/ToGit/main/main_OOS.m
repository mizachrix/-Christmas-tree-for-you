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

Input_param.eps = 0.4; % \epsilon
Input_param.Max_iter = 10000; % MaxIter
Input_param.tolerance = 1e-2; % \eta
Input_param.rho = 0.001; % \rho
Input_param.restrict = 0.02; % restriction at isolated approach

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

DeCC_TC = NaN(1,IR_sim);
CC_TC = NaN(1,IR_sim);
Gauss_TC = NaN(1,IR_sim);
RO_TC = NaN(1,IR_sim);

j = 1;

% For each coupled dataset, we pick N and N' samples
WPf_max = nWind(:,1:N_max,j)';
WPr_max = nWind(:,N_max+1:N_max+OOS_max,j)';
WPf = WPf_max(1:N,:);
WPr = WPr_max(1:OOS_sim,:);

    
% Build the corresponding data related to wind power production
all = [1:N];
system_info.Wscen = WPf(all,:)';
system_info.mu = mean(Wind)';%mean(wff)';%mean(system_info.Wscen,2); 
system_info.cov = cov(Wind);%cov(wff);   

temp2 = wff' - repmat(system_info.mu, 1, size(wff,1));
temp = temp2';
Nwind = size(system_info.Wmax,1);

% Loop for each individual run for 100 coupled datasets

for j = 1:IR_sim
    display('out of sample iteration:');
    j
    
    WPf_max = nWind(:,1:N_max,j)';
    WPr_max = nWind(:,N_max+1:N_max+OOS_max,j)';

    WPf = WPf_max(1:N,:);
    WPr = WPr_max(1:OOS_sim,:);
    
    % Build the corresponding data for Objective function
    all = [1:N];
    system_info.Wscen = WPf(all,:)';
    system_info.mu = mean(system_info.Wscen,2); 
    system_info.cov = cov(system_info.Wscen');
    
    % Build the corresponding data for RO
    system_info.Wscen_RO = de2bi(0:2^size(system_info.Wmax,1)-1)';
    system_info.Wexp_RO = mean(system_info.Wscen_RO,2); 
    system_info.xi = system_info.Wscen_RO - repmat(system_info.mu, 1, size(system_info.Wscen_RO,2));   
    system_info.xi2 = [ones(1,size(system_info.xi,2));system_info.xi];
    system_info.exp_xi_xit = [1,zeros(1,size(system_info.mu,1));zeros(size(system_info.mu,1),1),system_info.cov+zeros(size(system_info.mu,1),1)*zeros(1,size(system_info.mu,1))];
    system_info.mu2 = [1;zeros(size(system_info.mu,1),1)];%[1;si.mu];

    
    Coup_EL_Gas_CC = Coup_EL_Gas2_CC(system_info, Input_param);
    Coup_X_CC{j} = Coup_EL_Gas_CC.X;
    Coup_Y_CC{j} = Coup_EL_Gas_CC.Y;
    Coup_Q_CC{j} = Coup_EL_Gas_CC.Q;
    Coup_Obj_CC{j} = Coup_EL_Gas_CC.Obj;
    Coup_X_CC_real{j} = Coup_EL_Gas_CC.X(:,1) + Coup_EL_Gas_CC.X(:,2:Nwind+1) * system_info.Wscen;
    Coup_Y_CC_real{j} = Coup_EL_Gas_CC.Y(:,1) + Coup_EL_Gas_CC.Y(:,2:Nwind+1) * system_info.Wscen;
    Coup_Q_CC_real{j} = Coup_EL_Gas_CC.Q(:,1) + Coup_EL_Gas_CC.Q(:,2:Nwind+1) * system_info.Wscen;
    Coup_F_CC_real{j} = Coup_EL_Gas_CC.f + Coup_EL_Gas_CC.F * system_info.Wscen;
    Coup_CC_time(j) = toc;
    
    tic;
    
    Coup_EL_Gas_Gauss = Coup_EL_Gas2_Gauss(system_info, Input_param);
    Coup_X_Gauss{j} = Coup_EL_Gas_Gauss.X;
    Coup_Y_Gauss{j} = Coup_EL_Gas_Gauss.Y;
    Coup_Q_Gauss{j} = Coup_EL_Gas_Gauss.Q;
    Coup_Obj_Gauss{j} = Coup_EL_Gas_Gauss.Obj;
    Coup_X_Gauss_real{j} = Coup_EL_Gas_Gauss.X(:,1) + Coup_EL_Gas_Gauss.X(:,2:Nwind+1) * system_info.Wscen;
    Coup_Y_Gauss_real{j} = Coup_EL_Gas_Gauss.Y(:,1) + Coup_EL_Gas_Gauss.Y(:,2:Nwind+1) * system_info.Wscen;
    Coup_Q_Gauss_real{j} = Coup_EL_Gas_Gauss.Q(:,1) + Coup_EL_Gas_Gauss.Q(:,2:Nwind+1) * system_info.Wscen;
    Coup_F_Gauss_real{j} = Coup_EL_Gas_Gauss.f + Coup_EL_Gas_Gauss.F * system_info.Wscen;
    Coup_Gauss_time(j) = toc;
    
    tic;
    
    Coup_EL_Gas_Deter = Coup_EL_Gas2_Deter(system_info, Input_param);
    Coup_X_Deter{j} = Coup_EL_Gas_Deter.X;
    Coup_Y_Deter{j} = Coup_EL_Gas_Deter.Y;
    Coup_Q_Deter{j} = Coup_EL_Gas_Deter.Q;
    Coup_Obj_Deter{j} = Coup_EL_Gas_Deter.Obj;
    Coup_X_Deter_real{j} = Coup_EL_Gas_Deter.X(:,1) + Coup_EL_Gas_Deter.X(:,2:Nwind+1) * system_info.Wscen;
    Coup_Y_Deter_real{j} = Coup_EL_Gas_Deter.Y(:,1) + Coup_EL_Gas_Deter.Y(:,2:Nwind+1) * system_info.Wscen;
    Coup_Q_Deter_real{j} = Coup_EL_Gas_Deter.Q(:,1) + Coup_EL_Gas_Deter.Q(:,2:Nwind+1) * system_info.Wscen;
    Coup_F_Deter_real{j} = Coup_EL_Gas_Deter.f + Coup_EL_Gas_Deter.F * system_info.Wscen;
    Coup_Deter_time(j) = toc;
    
            
        tic;
        
        % Loop for each out-of-sample realization 
        for k = 1:OOS_sim
            system_info.Wreal = WPr(k,:)';
            system_info.DWreal = system_info.Wreal - system_info.mu;
                 

            % Solve real-time optimal power flow for the solution of CC
            % approximation
            RT_Coup_EL_Gas_CC = RT_solve_R(system_info,Coup_EL_Gas_CC.X(:,1),Coup_EL_Gas_CC.Y(:,1));
            CC_Obj_IR(j,k) = trace((RT_Coup_EL_Gas_CC.p_RT)'*system_info.A*(RT_Coup_EL_Gas_CC.p_RT)) + trace((RT_Coup_EL_Gas_CC.g_RT)'*system_info.Ag*(RT_Coup_EL_Gas_CC.g_RT)) + system_info.C'*(RT_Coup_EL_Gas_CC.p_RT) + system_info.Cg'*(RT_Coup_EL_Gas_CC.g_RT) + system_info.Clsh'*RT_Coup_EL_Gas_CC.lshed_RT;%RT_Coup_EL_Gas_CC.Obj_RT;  
            CC_lshed{j,k} = RT_Coup_EL_Gas_CC.lshed_RT;
            CC_x{j,k} = Coup_EL_Gas_CC.X(:,1) + Coup_EL_Gas_CC.X(:,2:Nwind+1) * system_info.DWreal;
            CC_y{j,k} = Coup_EL_Gas_CC.Y(:,1) + Coup_EL_Gas_CC.Y(:,2:Nwind+1) * system_info.DWreal;
            CC_q{j,k} = Coup_EL_Gas_CC.Q(:,1) + Coup_EL_Gas_CC.Q(:,2:Nwind+1) * system_info.DWreal;
            CC_f{j,k} =  Coup_EL_Gas_CC.f + Coup_EL_Gas_CC.F * system_info.DWreal;
            CC_flag(j,k) = RT_Coup_EL_Gas_CC.Flag;

            
            % Solve real-time optimal power flow for the solution of Gauss
            % approximation
            RT_Coup_EL_Gas_Gauss = RT_solve_R(system_info,Coup_EL_Gas_Gauss.X(:,1),Coup_EL_Gas_Gauss.Y(:,1));
            Gauss_Obj_IR(j,k) = trace((RT_Coup_EL_Gas_Gauss.p_RT)'*system_info.A*(RT_Coup_EL_Gas_Gauss.p_RT)) + trace((RT_Coup_EL_Gas_Gauss.g_RT)'*system_info.Ag*(RT_Coup_EL_Gas_Gauss.g_RT)) + system_info.C'*(RT_Coup_EL_Gas_Gauss.p_RT) + system_info.Cg'*(RT_Coup_EL_Gas_Gauss.g_RT) + system_info.Clsh'*RT_Coup_EL_Gas_Gauss.lshed_RT;%RT_Coup_EL_Gas_Gauss.Obj_RT;  
            Gauss_lshed{j,k} = RT_Coup_EL_Gas_Gauss.lshed_RT;
            Gauss_x{j,k} = Coup_EL_Gas_Gauss.X(:,1) + Coup_EL_Gas_Gauss.X(:,2:Nwind+1) * system_info.DWreal;
            Gauss_y{j,k} = Coup_EL_Gas_Gauss.Y(:,1) + Coup_EL_Gas_Gauss.Y(:,2:Nwind+1) * system_info.DWreal;
            Gauss_q{j,k} = Coup_EL_Gas_Gauss.Q(:,1) + Coup_EL_Gas_Gauss.Q(:,2:Nwind+1) * system_info.DWreal;
            Gauss_f{j,k} =  Coup_EL_Gas_Gauss.f + Coup_EL_Gas_Gauss.F * system_info.DWreal;
            Gauss_flag(j,k) = RT_Coup_EL_Gas_Gauss.Flag;
            
            % Solve real-time optimal power flow for the solution of Deter
            % approximation
            RT_Coup_EL_Gas_Deter = RT_solve_R(system_info,Coup_EL_Gas_Deter.X(:,1),Coup_EL_Gas_Deter.Y(:,1));
            Deter_Obj_IR(j,k) = trace((RT_Coup_EL_Gas_Deter.p_RT)'*system_info.A*(RT_Coup_EL_Gas_Deter.p_RT)) + trace((RT_Coup_EL_Gas_Deter.g_RT)'*system_info.Ag*(RT_Coup_EL_Gas_Deter.g_RT)) + system_info.C'*(RT_Coup_EL_Gas_Deter.p_RT) + system_info.Cg'*(RT_Coup_EL_Gas_Deter.g_RT) + system_info.Clsh'*RT_Coup_EL_Gas_Deter.lshed_RT;%RT_Coup_EL_Gas_Gauss.Obj_RT;  
            Deter_lshed{j,k} = RT_Coup_EL_Gas_Deter.lshed_RT;
            Deter_x{j,k} = Coup_EL_Gas_Deter.X(:,1) + Coup_EL_Gas_Deter.X(:,2:Nwind+1) * system_info.DWreal;
            Deter_y{j,k} = Coup_EL_Gas_Deter.Y(:,1) + Coup_EL_Gas_Deter.Y(:,2:Nwind+1) * system_info.DWreal;
            Deter_q{j,k} = Coup_EL_Gas_Deter.Q(:,1) + Coup_EL_Gas_Deter.Q(:,2:Nwind+1) * system_info.DWreal;
            Deter_f{j,k} =  Coup_EL_Gas_Deter.f + Coup_EL_Gas_Deter.F * system_info.DWreal;
            Deter_flag(j,k) = RT_Coup_EL_Gas_Deter.Flag;
                       

            
        end
        RT_time(j) = toc;
        
        % Calculation of expected cost
        
        if RT_Coup_EL_Gas_CC.Flag == 0
            CC_TC(j) = mean(CC_Obj_IR(j,:));
        else
            CC_TC(j) = NaN;
        end
        
        if RT_Coup_EL_Gas_Gauss.Flag == 0
            Gauss_TC(j) = mean(Gauss_Obj_IR(j,:));
        else
            Gauss_TC(j) = NaN;
        end
        
        if RT_Coup_EL_Gas_Deter.Flag == 0
            Deter_TC(j) = mean(Deter_Obj_IR(j,:));
        else
            Deter_TC(j) = NaN;
        end

    
end



Time = toc

CC_RC =mean(CC_TC)
Gauss_RC =mean(Gauss_TC)
Deter_RC =mean(Deter_TC)


Proj_CC = mean(cell2mat(Coup_Obj_CC))
Proj_Gauss = mean(cell2mat(Coup_Obj_Gauss))
Proj_Deter = mean(cell2mat(Coup_Obj_Deter))
