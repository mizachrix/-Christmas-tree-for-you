function[sol] = solve_X_ADMM_3_CC(si,input, Input_param)

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
    %Y = sdpvar(Ngas, Nwind+1, 'full'); % Linear decision rule for real-time gas production
    %Q = sdpvar(Npipes, Nwind+1, 'full'); % Linear decision rule for gas flow
    
%     s_obj = sdpvar(1, Nscen); % sigma variable for obj
%     
%     si.xi2 = [ones(1,size(si.xi,2));si.xi];
%     %si.exp_xi_xit = [1,si.mu';si.mu,si.cov+si.mu*si.mu'];
%     si.exp_xi_xit = [1,zeros(1,size(si.mu,1));zeros(size(si.mu,1),1),si.cov+zeros(size(si.mu,1),1)*zeros(1,size(si.mu,1))];
%     si.mu2 = [1;zeros(size(si.mu,1),1)];%[1;si.mu];
    CC_eps = Input_param.eps;
    CC_margin = round(sqrt((1-CC_eps)/(CC_eps)),2);
    
    % Constraints set
    CS = [];
    
    % Day-ahead constraints    
    %CS = [CS, si.Pmin <= p , p <= si.Pmax];
    %CS = [CS, si.Gmin <= g , g <= si.Gmax];
    
    CS = [CS, sum(X(:,1)) + sum(si.Wmax.*si.mu) - sum(si.D) == 0];
    CS = [CS, sum(X(:,2:Nwind+1), 1) == -si.Wmax'];
    %CS = [CS, Xc == diag(si.phi) * X];
    
%     CS = [CS, sum(g) - sum(si.PG*p) - sum(si.Dg) == 0];
%     CS = [CS, sum(Z, 1) == sum(si.PG * Y, 1)];
    %CS = [CS, si.GG * Y(:,1) - si.PG * X(:,1) - si.DG * si.Dg + si.PLG * Q(:,1) == 0];
    %CS = [CS, si.GG * Y(:,2:Nwind+1) == si.PG * X(:,2:Nwind+1) - si.PLG * Q(:,2:Nwind+1)];
    for i = 1:Nunits
        CS = [CS,  CC_margin * norm(X(i,2:Nwind+1) * si.cov^(1/2),2) <= (si.Pmax(i) - X(i,1)')];
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
    
%     CS = [CS, - repmat(si.FG, 1, Nscen) <= repmat(si.PLG * si.GG * g, 1, Nscen) - repmat(si.PLG * si.PG * p, 1, Nscen) -  repmat(si.PLG * si.DG * si.Dg, 1, Nscen) + si.PLG * si.GG * Z * si.xi - si.PLG * si.PG * Y * si.xi <= repmat(si.FG, 1, Nscen)];
%     CS = [CS, - repmat(si.FG, 1, Nscen) <=  Q * si.xi2 <= repmat(si.FG, 1, Nscen)];
    
    % Build the objective function 
    %Obj = 1/Nscen * sum(si.C'*X*si.xi2) + 1/Nscen * sum(si.Cg'*Y*si.xi2); 
    %Obj = trace(X'*si.A*X*si.exp_xi_xit) + trace(Y'*si.Ag*Y*si.exp_xi_xit) + si.C'*X*[1;si.mu] + si.Cg'*Y*[1;si.mu];
    Obj = trace(X'*si.A*X*si.exp_xi_xit) + si.C'*X*[1;zeros(size(si.mu,1),1)] - trace(input.lambda' * (input.K * X)) + input.rho/2 * trace( ((input.Xc) - (input.K * X))' * ((input.Xc) - (input.K * X)) );

    
    % Settings
    optim_options = sdpsettings('solver', 'gurobi','gurobi.TimeLimit',1000,'gurobi.NumericFocus',3,'verbose',0);

    % Solve
    sol = optimize(CS, Obj, optim_options);

    sol.X = value(X);
    %sol.Y = value(Y);
    %sol.Q = value(Q);
    sol.Obj = value(Obj);
    sol.Flag = sol.problem;
    
end