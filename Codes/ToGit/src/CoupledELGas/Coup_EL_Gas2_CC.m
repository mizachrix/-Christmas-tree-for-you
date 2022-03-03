function[sol] = Coup_EL_Gas2_CC(si, Input_param)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data-driven Distributed Operation of Electricity and Natural Gas Systems
    % Christos ORDOUDIS, Viet Anh NGUYEN, Jalal KAZEMPOUR, Pierre PINSON, Daniel KUHN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This script implements the Coupled Electricity and Natural gas model
    % with DR CC
    
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

    
    CC_eps = Input_param.eps;
    CC_margin = round(sqrt((1-CC_eps)/(CC_eps)),2);
    
    % Constraints set
    CS = [];
    
    % Day-ahead constraints    
    CS = [CS, sum(X(:,1)) + sum(si.Wmax.*si.mu) - sum(si.D) == 0];
    CS = [CS, sum(X(:,2:Nwind+1), 1) == -si.Wmax'];

    CS = [CS, sum(Y(:,1)) - sum(Xc(:,1)) - sum(si.Dg) == 0];
    CS = [CS, sum(Y(:,2:Nwind+1), 1) - sum(Xc(:,2:Nwind+1), 1) == 0];
    %CS = [CS,  X(:,2:Nwind+1) <= 0]; % Y(:,2:Nwind+1) <= 0,
    
    CS = [CS, si.GG * Y - si.PG * si.IG * Xc - si.DG * [si.Dg,zeros(size(si.Dg,1),Nwind)] + si.PLG * Q == 0];
    KK = diag(si.phi);
    K = KK(any(KK,2),:);
    CS = [CS, Xc == K * X];
    
    for i = 1:Nunits
        CS = [CS,  CC_margin * norm(X(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.Pmax(i) - X(i,1)')];
        CS = [CS,  CC_margin * norm(-X(i,2:Nwind+1) * si.cov^(1/2),2) <= (-si.Pmin(i) - (-X(i,1))')];
        CS = [CS, 0 <= si.Pmax(i) - X(i,1)'];
        CS = [CS, 0 <= -si.Pmin(i) - (-X(i,1))'];
    end
    
    for i = 1:Ngas
        CS = [CS,  CC_margin * norm(Y(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.Gmax(i) - Y(i,1)')];
        CS = [CS,  CC_margin * norm(-Y(i,2:Nwind+1) * si.cov^(1/2),2) <= (-si.Gmin(i) - (-Y(i,1))')];
        CS = [CS, 0 <= si.Gmax(i) - Y(i,1)'];
        CS = [CS, 0 <= -si.Gmin(i) - (-Y(i,1))'];        
    end
    
    for i = 1:Npipes
        CS = [CS,  CC_margin * norm(Q(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.FG(i) - Q(i,1)')];
        CS = [CS,  CC_margin * norm(-Q(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.FG(i) - (-Q(i,1))')];
        CS = [CS, 0 <= si.FG(i) - Q(i,1)'];
        CS = [CS, 0 <= si.FG(i) - (-Q(i,1))'];
    end
    
    for i = 1:Nlines
        CS = [CS,  CC_margin * norm(si.Qg(i,:) * X(:,2:Nwind+1) + si.Qw(i,:)*si.DiagWmax * si.cov^(1/2),2) <= (si.F(i) - (si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D)')];
        CS = [CS,  CC_margin * norm( -(si.Qg(i,:) * X(:,2:Nwind+1) + si.Qw(i,:)*si.DiagWmax) * si.cov^(1/2),2) <= (si.F(i) - (-(si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D))')];
        CS = [CS, 0 <= si.F(i) - (si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D)'];
        CS = [CS, 0 <= si.F(i) - (-(si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D))'];
    end
    

    % Build the objective function  
    Obj = trace(X'*si.A*X*si.exp_xi_xit) + trace(Y'*si.Ag*Y*si.exp_xi_xit) + si.C'*X*si.mu2 + si.Cg'*Y*si.mu2;

    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',2000,'gurobi.NumericFocus',3,'verbose',0,'gurobi.ScaleFlag',2,'gurobi.BarQCPConvTol',1e-1);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    sol.X = value(X);
    sol.Y = value(Y);
    sol.f = value(si.Qg * X(:,1) + si.Qw*si.DiagWmax*si.mu - si.Qd*si.D);
    sol.F = value(si.Qg * X(:,2:Nwind+1) + si.Qw*si.DiagWmax);
    sol.Q = value(Q);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;
    
end