function[sol] = DeCoup_Y_CC(si,input, Input_param)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data-driven Distributed Operation of Electricity and Natural Gas Systems
    % Christos ORDOUDIS, Viet Anh NGUYEN, Jalal KAZEMPOUR, Pierre PINSON, Daniel KUHN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This script solve the natural gas problem

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
    Xc = sdpvar(nnz(si.phi), Nwind+1, 'full'); % Linear decision rule for real-time power production 
    Y = sdpvar(Ngas, Nwind+1, 'full'); % Linear decision rule for real-time gas production
    Q = sdpvar(Npipes, Nwind+1, 'full'); % Linear decision rule for gas flow
    
    CC_eps = Input_param.eps;
    CC_margin = round(sqrt((1-CC_eps)/(CC_eps)),2);
   
    
    % Constraints set
    CS = [];
    
    % Day-ahead constraints    
   
    CS = [CS, si.GG * Y - si.PG * si.IG * input.K * input.X - si.DG * [si.Dg,zeros(size(si.Dg,1),Nwind)] + si.PLG * Q == 0];
    CS = [CS, sum(Y(:,1)) - sum(input.K * input.X(:,1)) - sum(si.Dg) == 0];
    CS = [CS, sum(Y(:,2:Nwind+1), 1) - sum(input.K * input.X(:,2:Nwind+1), 1) == 0];
    
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

    % Build the objective function 
    Obj = trace(Y'*si.Ag*Y*si.exp_xi_xit) + si.Cg'*Y*[1;zeros(size(si.mu,1),1)]; 
 
    
    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',1000,'gurobi.NumericFocus',3,'verbose',0,'gurobi.BarQCPConvTol',1e-1);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    %sol.X = value(X);
    sol.Y = value(Y);
    sol.Xc = value(Xc);
    sol.Q = value(Q);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;
    
end