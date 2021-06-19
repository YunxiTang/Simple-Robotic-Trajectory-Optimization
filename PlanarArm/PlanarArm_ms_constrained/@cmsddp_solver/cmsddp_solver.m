classdef cmsddp_solver < handle
    %MSDDPSOLVER constrained multiple shooting DDP solver
    
    properties
        M,             % shooting phase
        L,             % length of shooting phase
        Reg = 0.0,     % how much to regularize
        V_reg = 0.0,   % how much to regularize 
        Reg_Type = 2,  % 1->reg Quu (Default) / 2->reg Vxx
        eps = 1.0,     % eps: line-search parameter  
        gamma = 0.01,   % threshold to accept a FW step
        beta = 0.1,    % for line-search backtracking
        iter = 0,      % count iterations
        Jstore = []    % store real costs
        u_perturb = [],% constrol noise
        Lambda     ,   % dual variabels.      (path constraint) 
        Mu         ,   % penalty multipliers. (path constraint)
        Lambda_f,      % dual variabels.      (final constraint)
        Mu_f,          % penalty multipliers. (final constraint)
        phi = 2.0,      % penalty scaling parameter (path constraint) 
        phi_f = 2.5,    % penalty scaling parameter (final constraint)
        Path_Constraint,    % path constraints
        Final_Constraint,   % final constraints
        Cons_Vio = [],
        alpha_ = 1e-3
    end
    
    methods
        function obj = cmsddp_solver(path_constraint, final_constraint, params)
            % CMSDDP_SOLVER 
            disp('[INFO]: Calling Constrained Multiple Shooting DDP/SLQ solver.');
            obj.M = params.shooting_phase;
            obj.Reg_Type = params.Reg_Type;
            obj.L = params.N / obj.M + 1; 
            obj.Path_Constraint = path_constraint;
            obj.Final_Constraint = final_constraint;
            %%% initialize lambda and mu
            obj.Mu = cell(obj.M,1);
            obj.Lambda = cell(obj.M, 1);
            for k = 1:obj.M
                obj.Mu{k} = 2 * ones(path_constraint.n_ineq, obj.L-1);
                obj.Lambda{k} = 0 * ones(path_constraint.n_ineq, obj.L-1);
            end
            obj.Mu_f = 2 * ones(final_constraint.n_ineq, 1);
            obj.Lambda_f = 0 * ones(final_constraint.n_ineq, 1);
        end
        
        function [] = J_pushback(obj, J)
            % store costs
            obj.Jstore = [obj.Jstore J];
            
        end
        
        function [] = Update_iter(obj)
            % update iteration
            obj.iter = obj.iter + 1;
        end
        
        function [] = Update_Cons_Vio(obj,cons)
            obj.Cons_Vio = [obj.Cons_Vio cons];
        end
        
        function [] = dual_update(obj,xbar,ubar,params)
            % for path constraint
            for i=1:obj.M 
                for j=1:obj.L-1
                    x_ij = xbar{i}(:,j);
                    u_ij = ubar{i}(:,j);
                    obj.Lambda{i}(:,j) = max(0, obj.Lambda{i}(:,j) + obj.Mu{i}(:,j) .* obj.Path_Constraint.c(x_ij, u_ij));
                    obj.Mu{i}(:,j) = obj.phi .* obj.Mu{i}(:,j);
                end
            end
            
            % for final constraint
            x_end = xbar{obj.M}(:,end);
            obj.Lambda_f = obj.Lambda_f + obj.Mu_f .* obj.Final_Constraint.c(x_end);% max(0, obj.Lambda_f + obj.Mu_f .* obj.Final_Constraint.c(x_end));
            obj.Mu_f = obj.phi_f .* obj.Mu_f;
        end   
        
        function [] = solver_Callback(obj,xbar,ubar,params)
            % For plot
            nx = params.nx;
            dft = obj.CalDefect(xbar,params);
            figure(params.MapNo);
            for i = 1:obj.M
                clr = [i, 0.5 * i, 0.5 * i] / obj.M;
                figure(params.MapNo);
                
                plot(params.x0(1), params.x0(2), 'kp', 'MarkerFaceColor', 'b', 'MarkerSize', 15); hold on;
                plot(params.xf(1), params.xf(2), 'rh', 'MarkerFaceColor', 'r', 'MarkerSize', 15); hold on;
                
                plot(xbar{i}(1,:),xbar{i}(2,:),'Color',clr,'LineWidth',2.0);hold on;
                scatter(xbar{i}(1,1),xbar{i}(2,1),'MarkerFaceColor',clr); hold on;
                axis equal;
                grid on;
            end
            
            hold off;
            fig = figure(555);hold on;
            clf(fig);
            for k=1:nx
                title('Defects','Interpreter','latex','FontSize',20);
                subplot(nx,1,k);
                plot(dft(k,:),'s','LineWidth',2.0);
                hold off;
                grid on;
            end
        end
        
        function [J_idx,xsol,usol,max_cons] = simulate_phase(obj,rbt,cst,path_constraint,final_constraint,params,idx,x0,xbar,ubar,du,K)
            % simulate each inter-phase
            max_cons = 0;
            alpha = obj.eps;
            xsol = 0 * xbar;
            usol = 0 * ubar;
            J_idx = 0;
            xsol(:,1) = x0;
            xi = xsol(:,1);
            for i=1:(obj.L-1)
                dxi = xi - xbar(:,i);
                % Update with stepsize and feedback
                ui = ubar(:,i) + alpha*(du(:,i)) + K(:,:,i)*dxi;
                if params.clamp == 1
                    lb = params.umin * ones(params.nu, 1);
                    ub = params.umax * ones(params.nu, 1);
                    % clamp the control input
                    ui = max(lb, min(ub, ui));
                end
                usol(:,i) = ui;
                % Sum up costs
                lambda = obj.Lambda{idx}(:,i);
                Imu = obj.Path_Constraint.active_check(xi, ui, lambda, obj.Mu{idx}(:,i));
                J_idx = J_idx + cst.l_cost(xi, ui) + obj.Path_Constraint.AL_Term(xi, ui, lambda, Imu);
                % Propagate dynamics
                xi = rbt.rk(xi, ui, params.dt);
                % constraint violation
                c_i = path_constraint(xi, ui);
                if any(c_i > 0)
                    max_cons = max(max(c_i), max_cons);
                end
                xsol(:,i+1) = xi;
            end
            if idx==obj.M
                xN = xsol(:,end);
                lambdaN = obj.Lambda_f;
                ImuN = obj.Final_Constraint.active_check(xN, lambdaN, obj.Mu_f(:,end));
                J_idx = J_idx + cst.lf_cost(xsol(:,end)) + obj.Final_Constraint.AL_Term(xN, lambdaN, ImuN);
            end
            
        end
       
        function [J,xbar,ubar,Max_cons] = Init_Forward(obj,rbt,cst,path_constraint,final_constraint,params)
            % init forward simulation
            Max_cons = 0.0;
            xbar = cell(obj.M, 1);
            ubar = cell(obj.M, 1);
            J = 0;
            if isfield(params, {'xref'}) && params.warm_start == 1
                xbar = params.xref;
                ubar = params.uref;
            else
                for k = 1:obj.M
                    % make initial guess of intermediate states
                    xbar{k} = kron(ones(1, obj.L), params.xf);
                    ubar{k} = 0*ones(params.nu, obj.L-1);
                    xbar{k}(:,1) = params.x0 + (k-1) * (params.xf - params.x0) ./ obj.M;
                end
            end
            du = zeros(params.nu, obj.L-1);
            K  = -100*ones(params.nu, params.nx, obj.L-1);

            for i=1:obj.M
                x0 = xbar{i}(:,1);
                [J_idx,xbar{i},ubar{i},max_cons] = obj.simulate_phase(rbt,cst,path_constraint,final_constraint,params,i,x0,xbar{i},ubar{i},du,K);
                J = J + J_idx;
                Max_cons = max(Max_cons,max_cons);
            end
            obj.J_pushback(J);
        end
        
        
        function [dft] = CalDefect(obj,xbar,params)
            % compute all the defects
            dft = zeros(params.nx, params.shooting_phase-1);
            for i = 1:obj.M-1
                dft(:,i) = xbar{i}(:,end) - xbar{i+1}(:,1);
            end
        end
        
        function [J,x,u,Max_cons] = ForwardPass(obj,rbt,cst,path_constraint,final_constraint,params,xbar,ubar,du,K)
            % foward simulation
            
            Max_cons = 0;
            x = cell(obj.M, 1);
            u = cell(obj.M, 1);
            J = 0;
            [dft] = obj.CalDefect(xbar,params);
            new_dft = (1 - min(obj.eps,0.8)) .* dft;
            
            for k = 1:obj.M
                x{k} = 0 * xbar{k};
                u{k} = 0 * ubar{k};
            end
            for i=1:obj.M
                if i == 1
                    xbar{i}(:,1) = xbar{i}(:,1);
                    x0 = xbar{i}(:,1);
                else
                    %%% Method 1: From Crocoddyl
                    x0 = x{i-1}(:,end) -  new_dft(:,i-1);
                    
                    %%% Method 2: From Control Toolbox of ETHz
                    [fx, fu] = rbt.getLinSys(xbar{i-1}(:,end),ubar{i-1}(:,end),params.dt);
                    tilda = (fx + fu * K{i-1}(:,:,end)) * ((x{i-1}(:,end) - xbar{i-1}(:,end))) + fu * du{i-1}(:,end);
                    x0 = xbar{i}(:,1) + obj.eps * (tilda + dft(:,i-1));
                end
                [J_idx,x{i},u{i},max_cons] = obj.simulate_phase(rbt,cst,path_constraint,final_constraint,params,i,x0,xbar{i},ubar{i},du{i},K{i});
                J = J + J_idx;
                Max_cons = max(Max_cons,max_cons);
            end
            
        end
        
        function [dV,Vx,Vxx,du,K,success] = BackwardPass(obj,rbt,cst,xbar,ubar,dft,params)
            % backward propogation
            success = 1;
            % dV: tocompute the expected cost reduction
            dV = [0.0, 0.0];
            % initialize Vx, Vxx, du, K
            Vx = cell(obj.M, 1);
            Vxx = cell(obj.M, 1);
            du = cell(obj.M,1);
            K = cell(obj.M,1);
            for i = 1:obj.M
                du{i} = zeros(params.nu, obj.L-1);
                K{i} = zeros(params.nu, params.nx, obj.L-1);
                Vx{i} = zeros(params.nx, obj.L);
                Vxx{i} = zeros(params.nx, params.nx, obj.L);
            end
            xf = xbar{obj.M}(:,obj.L);
            [~,lfx,lfxx] = cst.lf_info(xf);
            
            %%%%%%%%%%% add constraints info%%%%%%%%%%%%%%%
            c_N = obj.Final_Constraint.c(xf);
            c_Nx = obj.Final_Constraint.algrad(xf);
            lambdaN = obj.Lambda_f;
            Imu_N = obj.Final_Constraint.active_check(xf, lambdaN, obj.Mu_f);
            Vx{obj.M}(:,obj.L) = lfx + c_Nx' * (lambdaN + Imu_N * c_N);
            Vxx{obj.M}(:,:,obj.L) = lfxx + c_Nx' * Imu_N * c_Nx;
            
            gap = zeros(params.nx,1);
            for i = (obj.M):-1:1
                if i < obj.M
                    gap = dft(:,i);
                    Vx{i}(:,obj.L) = Vx{i+1}(:,1) + Vxx{i+1}(:,:,1) * dft(:,i);
                    Vxx{i}(:,:,obj.L) = Vxx{i+1}(:,:,1);
                end
                for j = (obj.L-1):-1:1
                    if j < (obj.L-1)
                        gap = zeros(params.nx,1);
                    end
                    x_ij = xbar{i}(:,j);
                    u_ij = ubar{i}(:,j);
                    Imu_ij = obj.Path_Constraint.active_check(x_ij, u_ij, obj.Lambda{i}(:,j), obj.Mu{i}(:,j));
                    lambda_ij = obj.Lambda{i}(:,j);
                    Vx_ij = Vx{i}(:,j+1);
                    Vxx_ij = Vxx{i}(:,:,j+1);
                    [Qx,Qu,Qxx,Quu,Qux,Qxu,Quu_hat,Qux_hat] = cst.Qcms_info(rbt,cst,obj.Path_Constraint,lambda_ij,Imu_ij,x_ij,u_ij,Vx_ij,Vxx_ij,params,obj.V_reg);
                    % regularization
                    Quu_reg = Quu_hat + eye(params.nu)*obj.Reg;
                    % Make sure Quu is PD, if not, exit and increase regularization
                    [~, FLAG] = chol(Quu_reg-eye(params.nu)*1e-9);
                    if FLAG ~= 0 
                        % Quu is not PD, then increase Reg factor until Quu
                        % is PD
                        success = 0;
                        if params.Debug == 1
                            fprintf(' \t \t [SubSubInfo]: Reg=%5f:  Break Backward to Increase Reg. \n', obj.V_reg);
                        end
                        break
                    end
                    
                    % Standard Recursive Equations
                    if params.qp == 1 && obj.iter > 2
                        lb = params.umin * ones(params.nu, 1);
                        ub = params.umax * ones(params.nu, 1);
                        lower = lb - u_ij;
                        upper = ub - u_ij;
                        [kff,~,R,free] = boxQP(Quu_reg, Qu, lower, upper);
                        kfb = zeros(params.nu, params.nx);
                        if any(free)
                            Lfree = -R\(R'\Qux_hat(free,:));
                            kfb(free,:) = Lfree;
                        end
                    else
                        % without boxQP 
                        % more numerical stable
                        [R, ~] = chol(Quu_reg);
                        kff = -R\(R'\Qu);
                        kfb = -R\(R'\Qux_hat);
                    end
                    du{i}(:,j) = kff;
                    K{i}(:,:,j) = kfb;
                    Vx{i}(:,j)  = Qx + kfb'*Quu*kff + kfb'*Qu + Qxu*kff ;
                    Vxx{i}(:,:,j) = Qxx + Qxu*kfb + kfb'*Qux + kfb'*Quu*kfb;
                    delta1 = kff'*Qu + gap'*(Vx{i}(:,j) - Vxx{i}(:,:,j)*x_ij);
                    delta2 = kff'*Quu*kff + gap'*(2*Vxx{i}(:,:,j)*x_ij-Vxx{i}(:,:,j)*gap);
                    dV(1) = dV(1) + delta1;
                    dV(2) = dV(2) + 1/2 * delta2;
                end
                if success == 0
                    break
                end
            end
            if params.Debug == 1
                fprintf('\t \t [SubInfo]: Reg=%5f\n',obj.V_reg);
            end
        end
        
        function [V,x,u,max_cons,Flag] = ForwardIteration(obj,xbar,ubar,Vprev,consprev,du,K,DV,rbt,cst,path_constraint,final_constraint,params)
            % Check the Forward Pass is accepted or not (IMPORTANT!!!)
            % if not, adjust the line-search factor  (Armijo backtracking
            % line search)
            Flag = 1;
            V = 0;
            obj.eps = 1.0;
            alpha = obj.eps;
            while obj.alpha_ <= alpha
                % Try a step
                alpha = obj.eps;
                [V,x,u,max_cons] = obj.ForwardPass(rbt,cst,path_constraint,final_constraint,params,xbar,ubar,du,K);
                dV = alpha * DV(1) + alpha^2 * DV(2);
                ratio = (V - Vprev)/(dV);
                if params.Debug == 1
                    fprintf(' \t \t \t LS=%.3e \t A R=%.3e \t ER=%.3e \t Ratio=%.3e\n',...
                          alpha,V-Vprev, dV,ratio);
                end
                if dV < 0 && obj.gamma < ratio && ratio < 10 && max_cons <=  consprev
                    break
                elseif dV >= 0 && V < Vprev + 2 * dV && max_cons <=  consprev
                    break
                end
                % Else, do backtracking
                obj.eps = obj.beta * obj.eps;
            end
            if alpha < obj.alpha_
                V = Vprev;
                x = xbar;
                u = ubar;
                Flag = 0;
            end
            
        end
        
        function [xsol, usol, Ksol, dft, xbar, ubar] = Solve(obj,rbt,cst,path_constraint,final_constraint,params)
            % solve OCP
            % init rolling out
            [~,xbar,ubar,~] = obj.Init_Forward(rbt,cst,path_constraint,final_constraint,params);
            obj.Update_iter();
            % plot initial trajectories
            if params.plot == 1
                obj.solver_Callback(xbar,ubar,params);
            end
            [dft] = obj.CalDefect(xbar,params);

            [dV,~,~,du,K,~] = obj.BackwardPass(rbt,cst,xbar,ubar,dft,params);

            [Vbar,xbar,ubar,max_cons] = obj.ForwardPass(rbt,cst,path_constraint,final_constraint,params,xbar,ubar,du,K);
            obj.Update_iter();
            obj.Update_Cons_Vio(max_cons);
            [dft] = obj.CalDefect(xbar,params);
            consprev = max_cons;
            %%% start iteration
            while obj.iter < params.Max_iter
                if params.Debug == 1
                    fprintf('[INFO]: Iter. %3d   ||  Cost %.12e \n',obj.iter,Vbar);
                end
                % Set regularization back to 0 for next backward pass
                obj.V_reg = 1e-3;
                
                [dV,~,~,du,K,success] = obj.BackwardPass(rbt,cst,xbar,ubar,dft,params);
                
                Vprev = Vbar;

                %%% Forward Pass
                [Vbar,xbar,ubar,max_cons,Flag] = obj.ForwardIteration(xbar,ubar,Vprev,consprev,du,K,dV,rbt,cst,path_constraint,final_constraint,params);
                
                DELTA = 1;
                DELTA0 = 2;
                while (Flag == 0 || success == 0) 
                    DELTA = max(DELTA0, DELTA * DELTA0);
                    obj.V_reg = max(1e-3, obj.V_reg * DELTA);
                    if obj.V_reg > 1e3
                        break
                    end
                    [dV,~,~,du,K,success] = obj.BackwardPass(rbt,cst,xbar,ubar,dft,params);
                    [Vbar,xbar,ubar,max_cons,Flag] = obj.ForwardIteration(xbar,ubar,Vprev,consprev,du,K,dV,rbt,cst,path_constraint,final_constraint,params);
                end
                obj.Update_Cons_Vio(max_cons);
                % stop condition
                change = Vprev - Vbar;
                if abs(change) < params.stop  && obj.iter > 5
                    [xsol, usol, Ksol] = obj.assemble_solution(xbar, ubar, K, params);
                    fprintf('[INFO]: Exit with Convergence.\n');
                    return
                end
                consprev = max_cons;
                obj.Update_iter();
                [dft] = obj.CalDefect(xbar,params);
                if params.plot == 1 
                    obj.solver_Callback(xbar,ubar,params);
                end
                obj.J_pushback(Vbar);
                if mod(obj.iter,5) == 0
                    obj.dual_update(xbar,ubar,params);
                end
            end
            [xsol, usol, Ksol] = obj.assemble_solution(xbar, ubar, K, params);
            fprintf('[INFO]: Exit with Max Iterations.\n');
        end
        
        function [xsol, usol, Ksol] = assemble_solution(obj, xbar, ubar, K, params)
            %%% assemble the final solution
            xsol = zeros(params.nx, params.N+1);
            usol = zeros(params.nu, params.N);
            Ksol = zeros(params.nu, params.nx, params.N);
            if 1 < params.shooting_phase
                xsol(:,1:obj.L-1) = xbar{1}(:,1:(end-1));
                usol(:,1:obj.L-1) = ubar{1}(:,1:(end)); 
                Ksol(:,:,1:(obj.L-1)) = K{1}(:,:,1:(end));
                for k=2:params.shooting_phase                
                    if k < params.shooting_phase
                        xsol(:,(k-1)*(obj.L-1)+1:(k*(obj.L-1))) = xbar{k}(:,1:(end-1));
                    end
                    usol(:,(k-1)*(obj.L-1)+1:(k*(obj.L-1))) = ubar{k}(:,1:(end));
                    Ksol(:,:,(k-1)*(obj.L-1)+1:(k*(obj.L-1))) = K{k}(:,:,1:(end));
                end
                Me = params.shooting_phase;
                xsol(:,(Me-1)*(obj.L-1)+1:(Me*(obj.L-1)+1)) = xbar{k};
            else
                xsol = xbar{1};
                usol = ubar{1}(:,1:(end));
                Ksol(:,:,1:(obj.L-1)) = K{1}(:,:,1:(end));
            end  
        end
    end
end

