function[sol] = Coup_EL_Gas2_Deter(si, Input_param)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Energy and Reserves Dispatch with\\ Distributionally Robust Joint Chance Constraints
    % Christos ORDOUDIS, Viet Anh NGUYEN, Daniel KUHN, Pierre PINSON
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This script implements the Bonferroni approximation

    yalmip('clear')

    % Getting the number of thermals power plants, wind farms, scenarions,
    % transmission lines and nodes
    Nunits = size(si.Pmax,1);
    Nwind = size(si.Wmax,1);
    Nscen = size(si.Wscen_RO,2);
    Ngas = size(si.Gmax,1);
    Npipes = size(si.FG,1);
    Nlines = size(si.F,1);
    
    % Definition of variables
    X = sdpvar(Nunits, Nwind+1, 'full'); % Linear decision rule for real-time power production
    Xc = sdpvar(nnz(si.phi), Nwind+1, 'full'); % Linear decision rule for real-time power production 
    Y = sdpvar(Ngas, Nwind+1, 'full'); % Linear decision rule for real-time gas production
    Q = sdpvar(Npipes, Nwind+1, 'full'); % Linear decision rule for gas flow
%     tau1 = sdpvar(Nunits, 1);
%     tau2 = sdpvar(Nunits, 1);
%     
%     s_obj = sdpvar(1, Nscen); % sigma variable for obj
    
%     si.xi2 = [ones(1,size(si.xi,2));si.xi];
%     %si.exp_xi_xit = [1,si.mu';si.mu,si.cov+si.mu*si.mu'];
%     si.exp_xi_xit = [1,zeros(1,size(si.mu,1));zeros(size(si.mu,1),1),si.cov+zeros(size(si.mu,1),1)*zeros(1,size(si.mu,1))];
%     si.mu2 = [1;zeros(size(si.mu,1),1)];%[1;si.mu];
    CC_eps = Input_param.eps;
    CC_margin = 0;%round(norminv(1-CC_eps),2);
%     si.exp_xi_xit = round(si.exp_xi_xit,2);
%     si.mu2 = round(si.mu2,2);
%     si.mu = round(si.mu,2);
    
    % Constraints set
    CS = [];
    
    % Day-ahead constraints    
    %CS = [CS, si.Pmin <= p , p <= si.Pmax];
    %CS = [CS, si.Gmin <= g , g <= si.Gmax];
    
    CS = [CS, sum(X(:,1)) + sum(si.Wmax.*si.mu) - sum(si.D) == 0];
    CS = [CS, sum(X(:,2:Nwind+1), 1) == -si.Wmax'];

    CS = [CS, sum(Y(:,1)) - sum(Xc(:,1)) - sum(si.Dg) == 0];
    CS = [CS, sum(Y(:,2:Nwind+1), 1) - sum(Xc(:,2:Nwind+1), 1) == 0];
    CS = [CS,  X(:,2:Nwind+1) <= 0]; %  Y(:,2:Nwind+1) <= 0,
    
%     CS = [CS, sum(g) - sum(si.PG*p) - sum(si.Dg) == 0];
%     CS = [CS, sum(Z, 1) == sum(si.PG * Y, 1)];
%     CS = [CS, si.GG * Y(:,1) - si.PG * X(:,1) - si.DG * si.Dg + si.PLG * Q(:,1) == 0];
%     CS = [CS, si.GG * Y(:,2:Nwind+1) == si.PG * X(:,2:Nwind+1) - si.PLG * Q(:,2:Nwind+1)];
    
    CS = [CS, si.GG * Y - si.PG * si.IG * Xc - si.DG * [si.Dg,zeros(size(si.Dg,1),Nwind)] + si.PLG * Q == 0];
   % CS = [CS, Xc == diag(si.phi) * X];
    %CS = [CS, si.IG * Xc == diag(si.phi) * X];
    KK = diag(si.phi);
    K = KK(any(KK,2),:);
    CS = [CS, Xc == K * X];
    
%     for i = 1:Nunits
%         CS = [CS, sqrt(X(i,:) * si.exp_xi_xit * X(i,:)') <= CC_margin * (si.Pmax(i) - si.mu2'*X(i,:)')];
%         CS = [CS, sqrt(-X(i,:) * si.exp_xi_xit * (-X(i,:))') <= CC_margin * (si.Pmin(i) - si.mu2'*(-X(i,:))')];
%     end
    
%     for i = 1:Nunits
%         CS = [CS, CC_margin * (si.Pmax(i) - X(i,1)') == tau1(i)];
%         CS = [CS, 0 <= tau1(i)];
%         CS = [CS, norm(X(i,2:Nwind+1) * si.cov * X(i,2:Nwind+1)',2) <= tau1(i)];
%         
%         CS = [CS, CC_margin * (si.Pmin(i) - (-X(i,1))') == tau2(i)];
%         CS = [CS, 0 <= tau2(i)];
%         CS = [CS, norm(-X(i,2:Nwind+1) * si.cov * (-X(i,2:Nwind+1))',2) <= tau2(i)];        
%     end
   
    for i = 1:Nunits
        %CS = [CS,  CC_margin * norm(X(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.Pmax(i) - X(i,1)')];
        %CS = [CS,  CC_margin * norm(-X(i,2:Nwind+1) * si.cov^(1/2),2) <= (-si.Pmin(i) - (-X(i,1))')];
        CS = [CS, 0 <= si.Pmax(i) - X(i,1)'];
        CS = [CS, 0 <= -si.Pmin(i) - (-X(i,1))'];
    end
    
    for i = 1:Ngas
        %CS = [CS,  CC_margin * norm(Y(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.Gmax(i) - Y(i,1)')];
        %CS = [CS,  CC_margin * norm(-Y(i,2:Nwind+1) * si.cov^(1/2),2) <= (-si.Gmin(i) - (-Y(i,1))')];
        CS = [CS, 0 <= si.Gmax(i) - Y(i,1)'];
        CS = [CS, 0 <= -si.Gmin(i) - (-Y(i,1))'];        
    end
    
    for i = 1:Npipes
        %CS = [CS,  CC_margin * norm(Q(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.FG(i) - Q(i,1)')];
        %CS = [CS,  CC_margin * norm(-Q(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.FG(i) - (-Q(i,1))')];
        CS = [CS, 0 <= si.FG(i) - Q(i,1)'];
        CS = [CS, 0 <= si.FG(i) - (-Q(i,1))'];
    end
    
    for i = 1:Nlines
        %CS = [CS,  CC_margin * norm(si.Qg(i,:) * X(:,2:Nwind+1) + si.Qw(i,:)*si.DiagWmax * si.cov^(1/2),2) <= (si.F(i) - (si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D)')];
        %CS = [CS,  CC_margin * norm( -(si.Qg(i,:) * X(:,2:Nwind+1) + si.Qw(i,:)*si.DiagWmax) * si.cov^(1/2),2) <= (si.F(i) - (-(si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D))')];
        CS = [CS, 0 <= si.F(i) - (si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D)'];
        CS = [CS, 0 <= si.F(i) - (-(si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D))'];
    end
    
%     for i = 1:Nunits
%         Z1 = [1 + (X(i,2:Nwind+1) * si.cov * X(i,2:Nwind+1)'); 1 - (X(i,2:Nwind+1) * si.cov * X(i,2:Nwind+1)'); 2 * (CC_margin * (si.Pmax(i) - X(i,1)'))];
%         CS = [CS, cone(Z1)];
%         Z2 = [1 + (-X(i,2:Nwind+1) * si.cov * (-X(i,2:Nwind+1))'); 1 - (-X(i,2:Nwind+1) * si.cov * (-X(i,2:Nwind+1))'); 2 * (CC_margin * (si.Pmin(i) - (-X(i,1))'))];
%         CS = [CS, cone(Z2)];        
%     end
%     
%     for i = 1:Nunits
%         Z1 = [CC_margin * (si.Pmax(i) - si.mu2'*X(i,:)'); X(i,:) * si.exp_xi_xit^(1/2)];
%         CS = [CS, cone(Z1)];
%         Z2 = [CC_margin * (si.Pmin(i) - si.mu2'*(-X(i,:))'); -X(i,:) * si.exp_xi_xit^(1/2)];
%         CS = [CS, cone(Z2)];        
%     end
%     
%     CS = [CS, cone(Z)];
%     
%     for i = 1:Ngas
%         CS = [CS, sqrt(Y(i,:) * si.exp_xi_xit * Y(i,:)') <= CC_margin * (si.Gmax(i) - si.mu2'*Y(i,:)')];
%         CS = [CS, sqrt(-Y(i,:) * si.exp_xi_xit * (-Y(i,:))') <= CC_margin * (si.Gmin(i) - si.mu2'*(-Y(i,:))')];
%     end
%     
%     for i = 1:Nlines
%         CS = [CS, sqrt((si.Qg(i,:) * X + [si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D,si.Qw(i,:)*si.DiagWmax]) * si.exp_xi_xit * (si.Qg(i,:) * X + [si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D,si.Qw(i,:)*si.DiagWmax])') <= CC_margin * (si.F(i) - si.mu2'*(si.Qg(i,:) * X + [si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D,si.Qw(i,:)*si.DiagWmax])')];
%         CS = [CS, sqrt(-(si.Qg(i,:) * X + [si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D,si.Qw(i,:)*si.DiagWmax]) * si.exp_xi_xit * (-(si.Qg(i,:) * X + [si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D,si.Qw(i,:)*si.DiagWmax]))') <= CC_margin * (si.F(i) - si.mu2'* (-(si.Qg(i,:) * X + [si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D,si.Qw(i,:)*si.DiagWmax])'))];
%     end
        
   
%     CS = [CS, - repmat(si.F, 1, Nscen) <= si.Qg * X  * si.xi2 + repmat(si.Qw*si.DiagWmax*si.mu - si.Qd*si.D, 1, Nscen) + (si.Qw*si.DiagWmax) * si.xi <= repmat(si.F, 1, Nscen)];
    
%     CS = [CS, - repmat(si.FG, 1, Nscen) <= repmat(si.PLG * si.GG * g, 1, Nscen) - repmat(si.PLG * si.PG * p, 1, Nscen) -  repmat(si.PLG * si.DG * si.Dg, 1, Nscen) + si.PLG * si.GG * Z * si.xi - si.PLG * si.PG * Y * si.xi <= repmat(si.FG, 1, Nscen)];
    
%     for i = 1:Npipes
%         CS = [CS, sqrt(Q(i,:) * si.exp_xi_xit * Q(i,:)') <= CC_margin * (si.FG(i) - si.mu2'*Q(i,:)')];
%         CS = [CS, sqrt(-Q(i,:) * si.exp_xi_xit * (-Q(i,:))') <= CC_margin * (si.FG(i) - si.mu2'*(-Q(i,:))')];
%     end
    
    % Build the objective function 
    %Obj = 1/Nscen * sum(si.C'*X*si.xi2) + 1/Nscen * sum(si.Cg'*Y*si.xi2); 
    Obj = trace(X'*si.A*X*si.exp_xi_xit) + trace(Y'*si.Ag*Y*si.exp_xi_xit) + si.C'*X*si.mu2 + si.Cg'*Y*si.mu2;

    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',1000,'gurobi.NumericFocus',3,'verbose',0,'gurobi.ScaleFlag',2,'gurobi.BarQCPConvTol',1e-1);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    sol.X = value(X);
    sol.Y = value(Y);
    sol.Q = value(Q);
    sol.f = value(si.Qg * X(:,1) + si.Qw*si.DiagWmax*si.mu - si.Qd*si.D);
    sol.F = value(si.Qg * X(:,2:Nwind+1) + si.Qw*si.DiagWmax);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;
    
end