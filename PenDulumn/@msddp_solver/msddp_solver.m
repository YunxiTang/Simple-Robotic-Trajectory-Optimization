classdef msddp_solver < handle
    %MSDDPSOLVER multiple shooting DDP solver
    
    properties
        M,             % shooting phase
        L,             % length of shooting phase
        Reg = 0.0,     % how much to regularize
        Reg_Type = 1,  % 1->reg Quu (Default) / 2->reg Vxx
        eps = 1.0,     % eps: line-search parameter  
        gamma = 0.1,   % threshold to accept a FW step
        beta = 0.8,    % for line-search backtracking
        iter = 0,      % count iterations
        Jstore = []    % store real costs
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
        
        function [] = Update_iter(obj)
            % update iteration
            obj.iter = obj.iter + 1;
        end
        
        function [] = solver_Callback(obj,xbar,ubar,params)
            % For plot
            dft = obj.CalDefect(xbar,params);
            for i = 1:obj.M
                clr = [i, 0.5 * i, 0.5 * i] / obj.M;
                figure(111);
                plot(xbar{i}(1,:),xbar{i}(2,:),'Color',clr,'LineWidth',2.0);hold on;
                scatter(xbar{i}(1,1),xbar{i}(2,1),'MarkerFaceColor',clr); hold on;
            end
            xlabel('$\theta$','Interpreter','latex','FontSize',15);
            ylabel('$\dot \theta$','Interpreter','latex','FontSize',15);
            title('Phase Portrait','Interpreter','latex','FontSize',20);
            grid on;
            hold off;
            figure(555);hold on;
            title('Defects','Interpreter','latex','FontSize',20);
            subplot(2,1,1);
            plot(dft(1,:),'s','LineWidth',2.0);hold on;grid on;
            subplot(2,1,2);
            plot(dft(2,:),'s','LineWidth',2.0);hold on;grid on;
        end
        
        function [J_idx,xsol,usol] = simulate_phase(obj,rbt,cst,params,idx,x0,xbar,ubar,du,K)
            % simulate each inter-phase
            alpha = obj.eps;
            xsol = 0 * xbar;
            usol = 0 * ubar;
            J_idx = 0;
            xsol(:,1) = x0;
            xi = xsol(:,1);
            for i=1:(obj.L-1)
                dxi = xi - xbar(:,i);
                % Update with stepsize and feedback
                ui = ubar(:,i) + alpha*(du(:,i)+0/(1.5*obj.iter+1)*rand(1)) + K(:,:,i)*dxi;
                lb = params.umin * ones(params.nu, 1);
                ub = params.umax * ones(params.nu, 1);
                ui = max(lb, min(ub, ui));
                usol(:,i) = ui;
                % Sum up costs
                J_idx = J_idx + cst.l_cost(xi, ui);
                % Propagate dynamics
                xi = rbt.rk45(xi, ui, params.dt);
                xsol(:,i+1) = xi;
            end
            J_idx = J_idx + cst.lf_cost(xsol(:,end))*(idx==obj.M);
        end
       
        function [J,xbar,ubar] = Init_Forward(obj,rbt,cst,params)
            % init forward simulation
            xbar = cell(obj.M, 1);
            ubar = cell(obj.M, 1);
            J = 0;
            for k = 1:obj.M
                % make initial guess of intermediate states
                xbar{k} = zeros(params.nx, obj.L);
                ubar{k} = zeros(params.nu, obj.L);
                xbar{k}(:,1) = params.x0 + (k-1) * (params.xf - params.x0) / obj.M;
            end
            du = zeros(params.nu, obj.L);
            K  = zeros(params.nu, params.nx, obj.L);
            for i=1:obj.M
                x0 = xbar{i}(:,1);
                [J_idx,xbar{i},ubar{i}] = obj.simulate_phase(rbt,cst,params,i,x0,xbar{i},ubar{i},du,K);
                J = J + J_idx;
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
        
        function [J,x,u] = ForwardPass(obj,rbt,cst,params,xbar,ubar,du,K)
            % foward simulation
            x = cell(obj.M, 1);
            u = cell(obj.M, 1);
            J = 0;
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
                    x0 = x{i-1}(:,end) -  new_dft(:,i-1);
                    [fx_pre, fu_pre] = rbt.getLinSys(xbar{i-1}(:,end),ubar{i-1}(:,end));
                    [fx_n, fu_n] = rbt.getLinSys(xbar{i}(:,1),ubar{i}(:,1));
                    fx = (fx_pre * 2) / 2;
                    fu = (fu_pre * 2) / 2;
                    tilda = (fx + fu * K{i-1}(:,:,end)) * ((x{i-1}(:,end) - xbar{i-1}(:,end))) + fu * obj.eps*du{i-1}(:,end);
%                     x0 = xbar{i}(:,1) + 1.0*(tilda + dft(:,i-1));
                end
                [J_idx,x{i},u{i}] = obj.simulate_phase(rbt,cst,params,i,x0,xbar{i},ubar{i},du{i},K{i});
                J = J + J_idx;
            end
            
        end
        
        function [dV,Vx,Vxx,du,K,success] = BackwardPass(obj,rbt,cst,xbar,ubar,dft,params)
            % backward propogation
            success = 1;
            % dV: tocompute the expected cost reduction
            dV = 0.0;
            % initialize Vx, Vxx, du, K
            Vx = cell(obj.M, 1);
            Vxx = cell(obj.M, 1);
            du = cell(obj.M,1);
            K = cell(obj.M,1);
            for i = 1:obj.M
                du{i} = zeros(params.nu, obj.L);
                K{i} = zeros(params.nu, params.nx, obj.L);
                Vx{i} = zeros(params.nx, obj.L);
                Vxx{i} = zeros(params.nx, params.nx, obj.L);
            end
            xf = xbar{obj.M}(:,obj.L);
            [~,lfx,lfxx] = cst.lf_info(xf);
            Vx{obj.M}(:,obj.L) = lfx;
            Vxx{obj.M}(:,:,obj.L) = lfxx;
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
                    [Qx,Qu,Qxx,Quu,Qux,Qxu] = cst.Q_info(rbt,cst,x_ij,u_ij,Vx_ij,Vxx_ij,params);
                    % regularization
                    Quu_reg = Quu + eye(params.nu)*obj.Reg;
                    % Make sure Quu is PD, if not, exit and increase regularization
                    [~, FLAG] = chol(Quu_reg-eye(params.nu)*1e-9);
                    while FLAG ~= 0 
                        % Quu is not PD, then increase Reg factor until Quu
                        % is PD
                        obj.Reg = obj.Reg + 1e-3;
                        Quu_reg = Quu + eye(params.nu)*obj.Reg;
                        [~, FLAG] = chol(Quu_reg-eye(params.nu)*1e-9);
                        if params.Debug == 1
                            fprintf(' \t \t [SubSubInfo]: Reg=%5f:  Regularization <-- Increase Reg. \n', obj.Reg);
                        end
                    end
                    
                    % Standard Recursive Equations
                    kff = -Quu_reg\Qu;
                    kfb = -Quu_reg\Qux;
                    du{i}(:,j) = kff;
                    K{i}(:,:,j) = kfb;
                    Vx{i}(:,j)  = Qx + kfb'*Quu*kff + kfb'*Qu + Qxu*kff ;
                    Vxx{i}(:,:,j) = Qxx + Qxu*kfb + kfb'*Qux + kfb'*Quu*kfb;
                    alpha = obj.eps;
                    delta1 = kff'*Qu + gap'*(Vx{i}(:,j) - Vxx{i}(:,:,j)*x_ij);
                    delta2 = kff'*Quu*kff + gap'*(2*Vxx{i}(:,:,j)*x_ij-Vxx{i}(:,:,j)*gap);
%                     dV = dV + 1/2*kff'*Quu*kff + kff'*Qu;
                    dV = alpha * delta1 + 1/2*alpha^2*delta2;
                end
            end
            if params.Debug == 1
                fprintf('\t \t [SubInfo]: Reg=%5f\n',obj.Reg);
            end
        end
        
        function [V,x,u] = ForwardIteration(obj,xbar,ubar,Vprev,du,K,dV,rbt,cst,params)
            % Check the Forward Pass is accepted or not (IMPORTANT!!!)
            % if not, adjust the line-search factor  (Armijo backtracking
            % line search)
            V = 0;
            obj.eps = 1.0;
            alpha = obj.eps;
            while alpha > 1e-5
                % Try a step
                alpha = obj.eps;
                [V,x,u] = obj.ForwardPass(rbt,cst,params,xbar,ubar,du,K);
                ratio = (V - Vprev)/(dV);
                if params.Debug == 1
                    fprintf(' \t \t \t Alpha=%.3e \t Actual Reduction=%.3e \t Expected Reduction=%.3e \t Ratio=%.3e\n',...
                          alpha,V-Vprev, obj.gamma*dV,ratio);
                end
                if dV < 0 && V < Vprev + obj.gamma * dV
                    break
                elseif dV >= 0 && V < Vprev + 2 * dV
                    break
                end
                % Else, do backtracking
                obj.eps = obj.beta * obj.eps;
            end
            obj.J_pushback(V);
            if params.plot == 1 && mod(obj.iter,2) == 0
               obj.solver_Callback(x,u,params);
            end
        end
        
        function [xbar,ubar,du,K,dft] = Solve(obj,rbt,cst,params)
            % solve OCP
            % init rolling out
            [~,xbar,ubar] = obj.Init_Forward(rbt,cst,params);
            obj.Update_iter();
            % plot initial trajectories
            if params.plot == 1
                obj.solver_Callback(xbar,ubar,params);
            end
            [dft] = obj.CalDefect(xbar,params);

            [dV,Vx,Vxx,du,K,success] = obj.BackwardPass(rbt,cst,xbar,ubar,dft,params);
            
            %%% run a forward iteration
            if params.shooting_phase > 1
                % avoid closing gap at first iteration
                obj.eps = 0.9;
            end
            [Vbar,xbar,ubar] = obj.ForwardPass(rbt,cst,params,xbar,ubar,du,K);
            obj.Update_iter();
            [dft] = obj.CalDefect(xbar,params);
            
            %%% start iteration
            while obj.iter < params.Max_iter
                if params.Debug == 1
                    fprintf('[INFO]: Iteration %3d   ||  Cost %.12e \n',obj.iter,Vbar);
                end
                % run backward pass
                [dV,Vx,Vxx,du,K,success] = obj.BackwardPass(rbt,cst,xbar,ubar,dft,params);
                Vprev = Vbar;
                
                % Set regularization back to 0 for next backward pass
                obj.Reg = 0.0;
                %%% Forward Pass
                [Vbar,xbar,ubar] = obj.ForwardIteration(xbar,ubar,Vprev,du,K,dV,rbt,cst,params);
                obj.Update_iter();
                [dft] = obj.CalDefect(xbar,params);
                change = Vprev - Vbar;
                if change < params.stop 
                    break
                end
            end
        end
        function [xsol, usol, Ksol] = assemble_solution(obj, xbar, ubar, K, params)
            %%% assemble the final solution
            if 1 < params.shooting_phase
                xsol = xbar{1}(:,1:(end-1));
                usol = ubar{1}(:,1:(end-1));
                Ksol = zeros(params.nu, params.nx, params.N-1);
                Ksol(:,:,1:(obj.L-1)) = K{1}(:,:,1:(end-1));
                for k=2:params.shooting_phase                
                    if k < params.shooting_phase
                        xsol = [xsol xbar{k}(:,1:(end-1))];
                    end
                    usol = [usol ubar{k}(:,1:(end-1))];
                    Ksol(:,:,(k-1)*(obj.L-1)+1:(k*(obj.L-1))) = K{k}(:,:,1:(end-1));
                end
                xsol = [xsol xbar{k}];
            else
                xsol = xbar{1};
                usol = ubar{1}(:,1:(end-1));
                Ksol(:,:,1:(obj.L-1)) = K{1}(:,:,1:(end-1));
            end  
        end
    end
end

