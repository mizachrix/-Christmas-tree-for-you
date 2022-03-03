function [ sol ] = ADMM_CC( si, Input_param )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data-driven Distributed Operation of Electricity and Natural Gas Systems
    % Christos ORDOUDIS, Viet Anh NGUYEN, Jalal KAZEMPOUR, Pierre PINSON, Daniel KUHN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize X,Y and Lambda
    tic;
    
    input.lambda = zeros(size(si.IG,1),size(si.Wmax,1)+1);

    input.rho = Input_param.rho;
    Max_iter = Input_param.Max_iter;
    tolerance = Input_param.tolerance;
    
    KK = diag(si.phi);
    input.K = KK(any(KK,2),:);
    
    X_init = solve_X2_CC(si,input, Input_param);
    
    Y_init = solve_Y2_CC(si,input, Input_param);
    
    input.lambda = 1 + input.rho/2 * ( (Y_init.Xc) - (input.K * X_init.X) );
    input.X = X_init.X;
    input.Xc = Y_init.Xc;
    
    iter = 1;

    while iter <= Max_iter 

        % First, fix lambda, solve for X
        sol_X = solve_X_ADMM_3_CC(si,input, Input_param);
        
        input.X = sol_X.X;
        
        X_keep(:,:,iter) = sol_X.X;
        
        % Second, fix lambda, solve for Y
        sol_Y = solve_Y_ADMM_3_CC(si,input, Input_param);
        
        Y_keep(:,:,iter) = sol_Y.Y;
        Q_keep(:,:,iter) = sol_Y.Q;
        
        input.Xc = sol_Y.Xc;
        
        temp_lambda = input.lambda;
        
        % Update Lambda
        input.lambda = temp_lambda + input.rho * ( (sol_Y.Xc) - (input.K * sol_X.X) );
        

        dif_sub(:,:,iter) = (sol_Y.Xc) - (input.K * sol_X.X);
        
        Obj(iter) = trace(sol_X.X'*si.A*sol_X.X*si.exp_xi_xit) + trace(sol_Y.Y'*si.Ag*sol_Y.Y*si.exp_xi_xit) + si.C'*sol_X.X*[1;zeros(size(si.mu,1),1)] + si.Cg'*sol_Y.Y*[1;zeros(size(si.mu,1),1)];
        

     
        if iter > 1 && norm(squeeze(dif_sub(:,:,iter)),2) < tolerance && norm(input.rho*(input.K * X_keep(:,:,iter)-input.K * X_keep(:,:,iter-1)),2) < tolerance
            break
        end
        
        if iter > 1
            dual_res(:,:,iter-1) = input.rho*(input.K * X_keep(:,:,iter)-input.K * X_keep(:,:,iter-1));
        end
        
        iter = iter + 1;


    end
    Time = toc;
    
    %dif_plot = squeeze(dif_sub(:,5,:));
    
    sol.pr_res = dif_sub;
    sol.dual_res = dual_res;
    sol.X = X_keep;
    sol.Y = Y_keep;
    sol.Q = Q_keep;
    sol.Obj = Obj;
    sol.iter = iter;
    sol.Time = Time;
    

end

