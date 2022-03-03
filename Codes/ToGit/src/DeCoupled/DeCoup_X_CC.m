function[sol] = DeCoup_X_CC(si,input, Input_param)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data-driven Distributed Operation of Electricity and Natural Gas Systems
    % Christos ORDOUDIS, Viet Anh NGUYEN, Jalal KAZEMPOUR, Pierre PINSON, Daniel KUHN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This script solve the electricity problem

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
    
    CC_eps = Input_param.eps;
    CC_margin = round(sqrt((1-CC_eps)/(CC_eps)),2);
    
    % Constraints set
    CS = [];
    
    % Day-ahead constraints       
    CS = [CS, sum(X(:,1)) + sum(si.Wmax.*si.mu) - sum(si.D) == 0];
    CS = [CS, sum(X(:,2:Nwind+1), 1) == -si.Wmax'];

    for i = 1:Nunits
        CS = [CS,  CC_margin * norm(X(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.Pmax_dec(i) - X(i,1)')];
        CS = [CS,  CC_margin * norm(-X(i,2:Nwind+1) * si.cov^(1/2),2) <= (-si.Pmin(i) - (-X(i,1))')];
        CS = [CS, 0 <= si.Pmax(i) - X(i,1)'];
        CS = [CS, 0 <= -si.Pmin(i) - (-X(i,1))'];
    end  
    
    for i = 1:Nlines
        CS = [CS,  CC_margin * norm(si.Qg(i,:) * X(:,2:Nwind+1) + si.Qw(i,:)*si.DiagWmax * si.cov^(1/2),2) <= (si.F(i) - (si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D)')];
        CS = [CS,  CC_margin * norm( -(si.Qg(i,:) * X(:,2:Nwind+1) + si.Qw(i,:)*si.DiagWmax) * si.cov^(1/2),2) <= (si.F(i) - (-(si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D))')];
        CS = [CS, 0 <= si.F(i) - (si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D)'];
        CS = [CS, 0 <= si.F(i) - (-(si.Qg(i,:) * X(:,1) + si.Qw(i,:)*si.DiagWmax*si.mu - si.Qd(i,:)*si.D))'];
    end
    

    % Build the objective function 
    Obj = trace(X'*si.A_dec*X*si.exp_xi_xit) + si.C_dec'*X*[1;zeros(size(si.mu,1),1)];

    
    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',1000,'gurobi.NumericFocus',3,'verbose',0,'gurobi.BarQCPConvTol',1e-1);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    sol.X = value(X);
    %sol.Y = value(Y);
    %sol.Q = value(Q);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;
    
end