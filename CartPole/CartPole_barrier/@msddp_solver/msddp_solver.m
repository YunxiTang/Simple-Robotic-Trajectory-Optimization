classdef msddp_solver < handle
    %MSDDPSOLVER multiple shooting DDP solver (Barrier Function Version) 
    properties
        M,                   % shooting phase
        L,                   % length of shooting phase
        Reg = 0.0,          % how much to regularize
        V_reg = 1e-3,         % how much to regularize
        Reg_Type = 1,        % 1->reg Quu (Default) / 2->reg Vxx
        eps = 1.0,           % eps: line-search parameter  
        gamma = 0.01,        % threshold to accept a FW step
        beta = 0.1,          % for line-search backtracking
        iter = 1,            % count iterations
        i_LS = 0,            % count line search iter
        Jstore = [],         % store full costs
        J_real = [],         % store real costs
        Contract_Rate = []   % store constraction rate
        i_LSm = 10,          % max line search iterations 
        Cons_Vio = []
    end
    
    methods
        function obj = msddp_solver(params)
            % MSDDP_SOLVER 
            disp('[INFO]: Calling Multiple Shooting DDP/SLQ solver.');
            obj.M = params.shooting_phase;
            obj.Reg_Type = params.Reg_Type;
            obj.L = params.N / obj.M + 1; 
        end
        
        function [] = J_pushback(obj, J)
            % store costs
            obj.Jstore = [obj.Jstore J];
        end
        
        function [] = Jreal_pushback(obj, Jr)
            % store real costs
            obj.J_real = [obj.J_real Jr];
        end
        
        function [] = Rate_pushback(obj, rate)
            % store Contraction rate
            obj.Contract_Rate = [obj.Contract_Rate rate];
        end
        
        function [] = Update_iter(obj)
            % update iteration
            obj.iter = obj.iter + 1;
        end
        
        function [] = Update_Cons_Vio(obj,cons)
            obj.Cons_Vio = [obj.Cons_Vio cons];
        end
        
        function [] = solver_Callback(obj,xbar,ubar,params)
            % For plot
            %nx = params.nx;
            %dft = obj.CalDefect(xbar,params);
            for i = 1:obj.M
                clr = [i, 0.5 * i, 0.5 * i] / obj.M;
                figure(111);
                plot(xbar{i}(1,:),xbar{i}(2,:),'Color',clr,'LineWidth',2.0);
                hold on; 
            end
            axis equal;
            title('Phase Portrait','Interpreter','latex','FontSize',20);
            grid on;
            hold off;
            
            for i = 1:obj.M
                clr = [i, 0.5 * i, 0.5 * i] / obj.M;
                figure(222);
                plot(params.t{i}(1,1:end-1),ubar{i}(:,1:end),'Color',clr,'LineWidth',2.0);hold on; 
            end
            title('Control Input','Interpreter','latex','FontSize',20);
            grid on;
            hold off;
            
        end
        
        function [J_idx,Jr_idx,xsol,usol,max_cons] = simulate_phase(obj,rbt,cst,path_constraint,final_constraint,params,idx,x0,xbar,ubar,du,K)
            % simulate each inter-phase
            max_cons = 0;
            xref = params.xf;
            uref = zeros(params.nu, 1);
            alpha = obj.eps;
            xsol = 0 * xbar;
            usol = 0 * ubar;
            J_idx = 0;
            Jr_idx = 0;
            xsol(:,1) = x0;
            xi = xsol(:,1);
            for i=1:(obj.L-1)
                dxi = xi - xbar(:,i);
                % Update with stepsize and feedback
                ui = ubar(:,i) + alpha*(du(:,i)) + K(:,:,i)*dxi;
                
                % clamp the control input
                if params.clamp == 1
                    lb = params.umin * ones(params.nu, 1);
                    ub = params.umax * ones(params.nu, 1);
                    % clamp the control input
                    ui = max(lb, min(ub, ui));
                end
                usol(:,i) = ui;
                
                % Sum up costs
                J_idx = J_idx + cst.l_cost(xi, ui, xref, uref) + path_constraint.penalty(xi, ui)*params.dt;
                Jr_idx = Jr_idx + cst.l_cost(xi, ui, xref, uref);
                
                % Propagate dynamics
                xi = rbt.rk(xi, ui, params.dt);
                % constraint violation
                c_i = path_constraint.c(xi, ui);
                if any(c_i > 0)
                    max_cons = max(max(c_i), max_cons);
                end
                if norm(xi,2) > 1e3
                    % FP go unstable, stop the rollout.
                    % keep the x-u trajectories the same
                    % as last iteration
                    xsol = xbar;
                    usol = ubar;
                    J_idx = 1e10;
                    break
                end
                xsol(:,i+1) = xi;
            end
            Jr_idx = Jr_idx + (cst.lf_cost(xsol(:,end)))*(idx==obj.M);
            J_idx = J_idx + (cst.lf_cost(xsol(:,end))+final_constraint.penalty(xsol(:,end)))*(idx==obj.M);
        end
       
        function [J,Jr,xbar,ubar,Max_cons] = Init_Forward(obj,rbt,cst,path_constraint,final_constraint,params)
            % init forward simulation
            Max_cons = 0.0;
            xbar = cell(obj.M, 1);
            ubar = cell(obj.M, 1);
            J = 0;
            Jr = 0;
            if isfield(params, {'x_warm','u_warm'})
                % utilize the warm-start trajectories
                xbar = params.x_warm;
                ubar = params.u_warm;
            else
                for k = 1:obj.M
                    % make initial guess of intermediate states
                    xbar{k} = kron(ones(1, obj.L), params.xf);
                    ubar{k} = -ones(params.nu, obj.L-1);
                    xbar{k}(:,1) = params.x0 + (k-1) * (params.xf - params.x0) ./ obj.M;
                end
            end
            du = zeros(params.nu, obj.L-1);
            K  = zeros(params.nu, params.nx, obj.L-1);
            for i=1:obj.M
                x0 = xbar{i}(:,1);
                [J_idx,Jr_idx,xbar{i},ubar{i},max_cons] = obj.simulate_phase(rbt,cst,path_constraint,final_constraint,params,i,x0,xbar{i},ubar{i},du,K);
                J = J + J_idx;
                Jr = Jr + Jr_idx;
                Max_cons = max(Max_cons,max_cons);
            end
            obj.J_pushback(J);
            obj.Jreal_pushback(Jr);
        end
        
        
        function [dft] = CalDefect(obj,xbar,params)
            % compute the defects
            dft = zeros(params.nx, params.shooting_phase-1);
            for i = 1:obj.M-1
                dft(:,i) = xbar{i}(:,end) - xbar{i+1}(:,1);
            end
        end
        
        function [J,Jr,x,u,Max_cons] = ForwardPass(obj,rbt,cst,path_constraint,final_constraint,params,xbar,ubar,du,K)
            % foward simulation
            Max_cons = 0;
            x = cell(obj.M, 1);
            u = cell(obj.M, 1);
            J = 0;
            Jr = 0;
            [dft] = obj.CalDefect(xbar,params);
            new_dft = (1 - obj.eps) .* dft;
            
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
%                     [fx, fu] = rbt.getLinSys(xbar{i-1}(:,end),ubar{i-1}(:,end), params.dt);
%                     tilda = (fx + fu * K{i-1}(:,:,end)) * ((x{i-1}(:,end) - xbar{i-1}(:,end))) + fu * obj.eps * du{i-1}(:,end);
%                     x0 = xbar{i}(:,1) +  (1.0) * (tilda) + 1.0 * dft(:,i-1);
                end
                [J_idx,Jr_idx,x{i},u{i},max_cons] = obj.simulate_phase(rbt,cst,path_constraint,final_constraint,params,i,x0,xbar{i},ubar{i},du{i},K{i});
                J = J + J_idx;
                Jr = Jr + Jr_idx;
                Max_cons = max(Max_cons,max_cons);
            end
            
        end
        
        function [dV,Vx,Vxx,du,K,Grad,success] = BackwardPass(obj,rbt,cst,path_constraint,final_constraint,xbar,ubar,dft,params)
            % backward propogation
            success = 1;
            Grad  = [];
            xref = params.xf;
            uref = zeros(params.nu, 1);
            % dV: to compute the expected cost reduction
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
            %%%%
            [~,lfx,lfxx] = cst.lf_info(xf);
            [~, Penal_x, Penal_xx] = final_constraint.penalty_info(xf);
            Vx{obj.M}(:,obj.L) = lfx + Penal_x;
            Vxx{obj.M}(:,:,obj.L) = lfxx + Penal_xx;
            %%%%
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
                    Vx_ij = Vx{i}(:,j+1);
                    Vxx_ij = Vxx{i}(:,:,j+1);
                    [Qx,Qu,Qxx,Quu,Qux,Qxu,Quu_hat,Qux_hat] = cst.Qcs_info(rbt,cst,x_ij,u_ij,xref,uref,Vx_ij,Vxx_ij,params,path_constraint,obj.V_reg);
                    Grad = [Grad norm(Qu,2)];
                    % regularization
                    Quu_reg = Quu_hat + eye(params.nu)*obj.Reg;
                    % Make sure Quu is PD, if not, exit and increase regularization
                    [~, FLAG] = chol(Quu_reg-eye(params.nu)*1e-6);
                    if FLAG > 0 
                        % Quu is not PD, then increase Reg factor until Quu
                        % is PD
                        success = 0;
                        if params.Debug == 1
                            fprintf(' \t \t [SubSubInfo]: Non PD Qxx! Reg=%5f.\n', obj.V_reg);
                        end
                        return
                    end
                    
                    % Standard Recursive Equations
                    % add boxQP solver here
                    if params.qp == 1 && obj.iter > 5
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
                        % use cholesky decomposition, more numerical stable
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
            end
            if params.Debug == 1
                fprintf('\t \t [SubInfo]: Reg=%5f\n',obj.V_reg);
            end
        end
        
        function [V,Vr,x,u,max_cons,Flag] = ForwardIteration(obj,xbar,ubar,Vprev,consprev,du,K,DV,rbt,cst,path_constraint,final_constraint,params)
            % Check the Forward Pass is accepted or not (IMPORTANT!!!)
            % if not, adjust the line-search factor  (Armijo backtracking
            % line search)
            V = 0;
            Vr = 0;
            obj.eps = 1.0;
            Flag = 1;
            alpha = obj.eps;
            while alpha > 1e-3
                % Try a step
                alpha = obj.eps;
                [V,Vr,x,u,max_cons] = obj.ForwardPass(rbt,cst,path_constraint,final_constraint,params,xbar,ubar,du,K);
                dV = alpha * DV(1) + alpha^2 * DV(2);
                ratio = (V - Vprev)/(dV);
                if params.Debug == 1
                    fprintf(' \t \t \t LS=%.3e \t AR=%.3e \t ER=%.3e \t Ratio=%.3e\n',...
                              alpha,V-Vprev, dV,ratio);
                end
                if obj.iter == 0
                    break
                end
                if dV < 0 && V < Vprev + obj.gamma * dV && max_cons <= 1.2* consprev
                    break
                elseif dV >= 0 && V < Vprev + 2 * dV && max_cons <= 1.2* consprev
                    break
                end
                % Else, do backtracking
                obj.eps = obj.beta * obj.eps;
            end
            if alpha <= 1e-3
                V = Vprev;
                x = xbar;
                u = ubar;
                Flag = 0;
            end
        end
        
        function [xsol, usol, Ksol, xbar, ubar] = Solve(obj,rbt,cst,path_constraint,final_constraint,params)
            % solve OCP
            % init rolling out
            [Vbar,~,xbar,ubar,max_cons] = obj.Init_Forward(rbt,cst,path_constraint,final_constraint,params);
            consprev = max_cons;
            % plot initial trajectories
            if params.plot == 1
                obj.solver_Callback(xbar,ubar,params);
            end
            [dft] = obj.CalDefect(xbar,params);
            
            while obj.iter <= params.Max_iter
                if params.Debug == 1
                    fprintf('[INFO]: Iter. %3d   ||  Cost %.12e ||  Cons. Vio. %.12e\n',obj.iter, Vbar, max_cons);
                end
                
                % Set regularization back to 0 for next backward pass
                obj.V_reg = 1e-6;
                obj.Reg = 1 / (exp(2*obj.iter)-2);
                % run backward pass
                [dV,~,~,du,K,Grad,success] = obj.BackwardPass(rbt,cst,path_constraint,final_constraint,xbar,ubar,dft,params);
                
                Vprev = Vbar;
                
                %%% Forward Pass
                [Vbar,Vr,xbar,ubar,max_cons,Flag] = obj.ForwardIteration(xbar,ubar,Vprev,consprev,du,K,dV,rbt,cst,path_constraint,final_constraint,params);

                DELTA = 5;
                DELTA0 = 5;
                
                obj.i_LS = 0;
                while (Flag == 0 || success == 0) && obj.i_LS <= obj.i_LSm
                    DELTA = max(DELTA0, DELTA * DELTA0);
%                     obj.V_reg = max(1e-6, obj.V_reg * DELTA);
                    obj.V_reg = 100 * obj.V_reg;
                    obj.i_LS = obj.i_LS + 1;
                    if obj.V_reg > 1e3
                          path_constraint.update_t();
                          final_constraint.update_t();
                          if params.Debug == 1
                            disp('~~~~~~~~~~~~~~~~LS & Reg Fialed. Update RLB.~~~~~~~~~~~~~~~~~~~');
                          end
                          [Vbar,Vr,xbar,ubar,max_cons] = obj.ForwardPass(rbt,cst,path_constraint,final_constraint,params,xbar,ubar,du,K);
                          Vprev = Vbar;
                          obj.V_reg = 1e-6;
                          DELTA = 5;
                          continue
                    end
                    [dV,~,~,du,K,Grad,success] = obj.BackwardPass(rbt,cst,path_constraint,final_constraint,xbar,ubar,dft,params);
                    [Vbar,Vr,xbar,ubar,max_cons,Flag] = obj.ForwardIteration(xbar,ubar,Vprev,consprev,du,K,dV,rbt,cst,path_constraint,final_constraint,params); 
                end
                
                if params.plot == 1 && mod(obj.iter,1) == 0
                    obj.solver_Callback(xbar,ubar,params);
                end
                [dft] = obj.CalDefect(xbar,params);
                V_change = Vprev - Vbar;
                DU = cell2mat(du);
                DU = reshape(DU,(params.nu*params.shooting_phase*(obj.L-1)),1);
                
                % Check the exit criteria
                if obj.iter > 2
                    if ((V_change) < params.stop) && max_cons <= 1e-6
                        fprintf('[INFO]: Value Function Converge!');
                        break
                    elseif all(DU<=1e-3) && max_cons <= 1e-6
                        fprintf('[INFO]: Policy Function Converge!');
                        break
                    elseif max(Grad)<=1e-8  && max_cons <= 1e-6
                        fprintf('[INFO]: Gradient Converge!');
                        break
                    elseif consprev - max_cons < 0.0
                        fprintf('[INFO]: Const. Satisfication!');
                        break
                    end
                end
                obj.J_pushback(Vbar);
                obj.Jreal_pushback(Vr);
                obj.Update_iter();
                obj.Update_Cons_Vio(max_cons);
                consprev = max_cons;
                % Update RLB functions
                if  (V_change) < 0.1 %mod(obj.iter, 10)==0
                    path_constraint.update_t();
                    final_constraint.update_t();
                    disp('~~~~~~~~~~~~~~~~Update RLB.~~~~~~~~~~~~~~~~~~~');
                end
            end
            [xsol, usol, Ksol] = obj.assemble_solution(xbar, ubar, K, params);
        end
        
        function [xsol, usol, Ksol] = assemble_solution(obj, xbar, ubar, K, params)
            %%% Assemble the final solution
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

