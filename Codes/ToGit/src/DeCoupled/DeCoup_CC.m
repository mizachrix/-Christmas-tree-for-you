function [ sol ] = DeCoup_CC( si, Input_param )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data-driven Distributed Operation of Electricity and Natural Gas Systems
    % Christos ORDOUDIS, Viet Anh NGUYEN, Jalal KAZEMPOUR, Pierre PINSON, Daniel KUHN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize alpha and beta
    tic;
    
    
    Max_iter = Input_param.Max_iter;
    
    restrict = Input_param.restrict;
    restrict_v = ones(1,size(si.phi2,1)) * restrict;
    logic = [1 0 0 0 1 0 1 0 0 0 1 0];
    restrict_l = zeros(1,size(si.phi,2));
    restrict_l(logic == 1) = restrict_v;
    restrict_l = restrict_l';
    Nwind = size(si.Wmax,1);
    KK = diag(si.phi);
    input.K = KK(any(KK,2),:);
    
    si.C_dec = (si.phi * mean(si.Cg))' + si.C;
    si.A_dec = diag(0.0001*si.C_dec);
    si.Pmax_dec = si.Pmax;
    
    iter = 1;

    while iter <= Max_iter 

        % First, fix lambda, solve for X
        sol_X = DeCoup_X_CC(si,input, Input_param);
        
        input.X = sol_X.X;
        
        X_keep(:,:,iter) = sol_X.X;
        
        % Second, fix lambda, solve for Y
        sol_Y = DeCoup_Y_CC(si,input, Input_param);
        
        Y_keep(:,:,iter) = sol_Y.Y;
        Q_keep(:,:,iter) = sol_Y.Q;
        
        Obj(iter) = trace(sol_X.X'*si.A*sol_X.X*si.exp_xi_xit) + trace(sol_Y.Y'*si.Ag*sol_Y.Y*si.exp_xi_xit) + si.C'*sol_X.X*[1;zeros(size(si.mu,1),1)] + si.Cg'*sol_Y.Y*[1;zeros(size(si.mu,1),1)];
        
        if sol_Y.Flag == 0
            break
        end
        
        si.Pmax_dec = si.Pmax_dec - si.Pmax .* restrict_l;
          
        iter = iter + 1;
        

    end
    
    Time = toc;
    
    sol.X = sol_X.X;
    sol.Y = sol_Y.Y;
    sol.Q = sol_Y.Q;
    sol.f = value(si.Qg * sol_X.X(:,1) + si.Qw*si.DiagWmax*si.mu - si.Qd*si.D);
    sol.F = value(si.Qg * sol_X.X(:,2:Nwind+1) + si.Qw*si.DiagWmax);
    sol.Obj = Obj;
    sol.iter = iter;
    sol.Time = Time;
    
    
end

