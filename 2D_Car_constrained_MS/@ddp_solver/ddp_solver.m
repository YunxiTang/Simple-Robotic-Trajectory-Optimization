
classdef ddp_solver < handle
    %DDPSOLVER ddp solver
    
    properties
        Reg = 0.000,   % how much to regularize
        Reg_Type = 1,  % 1->reg Quu (Default) / 2->reg Vxx
        eps = 1.0,     % eps: line-search parameter  
        gamma = 0.1,   % threshold to accept a FW step
        beta = 0.8,    % for line-search backtracking
        iter = 0,      % count iterations
        Jstore = [],   % store real costs
        u_perturb = []
    end
    
    methods
        function obj = ddp_solver(params)
            %DDP_SOLVER 
            disp('[INFO]: Calling DDP/SLQ solver.');
            if nargin > 0
                obj.Reg_Type = params.Reg_Type;
                rng('default');
                obj.u_perturb = rand(params.nu, params.N-1);
            end
        end
        
        function [] = J_pushback(obj, J)
            obj.Jstore = [obj.Jstore J];
        end
        
        function [] = Update_iter(obj)
            obj.iter = obj.iter + 1;
        end
        
        function [] = solver_Callback(obj,xbar,ubar,params)
            %%% for plot
            fig = figure(111);
            clf(fig);
            nx = params.nx;
            state_clr = ['r','g','b','k'];
            for i=1:nx
                plot(params.tax,xbar(i,:),'Color',state_clr(i),'LineWidth',2.0);hold on;
                grid on;
            end
            legend("$x$","$y$","$\theta$","$v$",'Interpreter','latex','FontSize',12);
            
            figure(222);
            plot(xbar(1,:),xbar(2,:),'k-','LineWidth',1.0);hold on;
            plot(xbar(1,1:5:end),xbar(2,1:5:end),'ro','LineWidth',2.0,'MarkerSize',8.0);hold off;
            grid on;
        end
        
        function [xbar, ubar] = Init_Forward(obj,rbt,params)
            %METHOD1 Init rollout
            ubar = zeros(params.nu, params.N);
            xbar = zeros(params.nx, params.N);

            xbar(:,1) = params.x0;
            xi = params.x0;
            % Make initial Guess of Trajectory
            for i = 1:(params.N-1)
                % Option 1: PD Control
%                 ui = -[15 5] * (xi - params.xf);
                % Option 2: Zero Control
                ui = zeros(params.nu,1);
                % Option 3: Random Control
                ubar(:,i) = ui;
                xi = rbt.rk(xi,ui,params.dt);
                xbar(:,i+1) = xi;
            end
            obj.solver_Callback(xbar,ubar,params);
        end
        
        function [J,x,u] = ForwardPass(obj,rbt,cst,params,xbar,ubar,du,K)
            %%% forward rollout
            alpha = obj.eps;
            J = 0;
            x = 0*xbar;
            u = 0*ubar;
            x(:,1) = xbar(:,1);
            xi = xbar(:,1);
            for i=1:params.N-1
                dxi = xi - xbar(:,i);
                % Update with stepsize and feedback
                % add noise to action-space (faster convergence for simple system)
                ui = ubar(:,i) + alpha*(du(:,i)+ 0.0 / obj.iter * obj.u_perturb(:,i)) + K(:,:,i)*dxi;
%                 if abs(ui) > params.umax
%                     ui = sign(ui)*params.umax;
%                 end
                u(:,i) = ui;
                % Sum up costs
                J = J + cst.l_cost(xi, ui);
                % Propagate dynamics
                xi = rbt.rk(xi, ui, params.dt);
                x(:,i+1) = xi;
            end
            J = J + cst.lf_cost(x(:,params.N));
            if params.plot == 1 && mod(obj.iter,1) == 0
                obj.solver_Callback(xbar,ubar,params);
            end
        end
        
        function [dV,Vx,Vxx,du,K,success] = BackwardPass(obj,rbt,cst,xbar,ubar,params)
            %%% Backward pass
            success = 1;
            % Initialization
            % compute the expected cost reduction
            dV = 0.0;
            Vx = zeros(params.nx, params.N);
            Vxx = zeros(params.nx, params.nx, params.N);
            du = zeros(params.nu, params.N);
            K = zeros(params.nu, params.nx, params.N);
            
            xf = xbar(:,params.N);
            [~,lfx,lfxx] = cst.lf_info(xf);
            Vx(:,params.N) = lfx;
            Vxx(:,:,params.N) = lfxx;
            for i = (params.N-1):-1:1
                xi = xbar(:,i);
                ui = ubar(:,i);
                Vxi = Vx(:,i+1);
                Vxxi = Vxx(:,:,i+1);
                [Qx,Qu,Qxx,Quu,Qux,Qxu] = cst.Q_info(rbt,cst,xi,ui,Vxi,Vxxi,params);
                % regularization
                Quu_reg = Quu + eye(params.nu)*obj.Reg;
                % Make sure Quu is PD, if not, exit and increase regularization
                [~, FLAG] = chol(Quu_reg-eye(params.nu)*1e-9);
                if FLAG ~= 0 
                    % Quu is not PD, then break out to increase Reg factor
                    if params.Debug == 1
                        fprintf(' \t \t [SubSubInfo]: Break BackWardPass And Increase Reg. \n');
                    end
                    success = 0;
                    break
                end
                % Standard Recursive Equations
                kff = -Quu_reg\Qu;
                kfb = -Quu_reg\Qux;
                du(:,i)  = kff;
                K(:,:,i) = kfb;
                Vx(:,i)  = Qx + kfb'*Quu*kff + kfb'*Qu + Qxu*kff ;
                Vxx(:,:,i) = Qxx + Qxu*kfb + kfb'*Qux + kfb'*Quu*kfb;
                dV = dV + 1/2*kff'*Quu*kff + kff'*Qu;
            end   
        end
        
        function [V,x,u] = ForwardIteration(obj,xbar,ubar,Vprev,du,K,dV,rbt,cst,params)
            % Check the Forward Pass is accepted or not (IMPORTANT!!!)
            % if not, adjust the line-search factor  (Armijo backtracking
            % line search)
            V = 0;
            obj.eps = 1.0;
            alpha = obj.eps;
            x = 0*xbar;
            u = 0*ubar;
            while alpha > 1e-5
                % Try a step
                alpha = obj.eps;
                [V,x,u] = obj.ForwardPass(rbt,cst,params,xbar,ubar,du,K);
                ratio = (V - Vprev)/(alpha*(1-alpha/2)*dV);
                if params.Debug == 1
                    fprintf(' \t \t \t Alpha=%.3e \t Actual Reduction=%.3e \t Expected Reduction=%.3e \t Ratio=%.3e\n',...
                          alpha,V-Vprev, obj.gamma*alpha*(1-alpha/2)*dV,ratio);
                end
                if V < Vprev + obj.gamma * alpha*(1-alpha/2)*dV
                    break
                end
                % Else, do backtracking
                obj.eps = obj.beta * obj.eps;
            end
        end
        
        function [xbar, ubar, K, du] = Solve(obj,rbt,cst,params)
            % init rolling out
            [xbar, ubar] = obj.Init_Forward(rbt,params);
            obj.Update_iter();
            du = zeros(params.nu, params.N);
            K = zeros(params.nu, params.nx, params.N);
            [Vbar,xbar,ubar] = obj.ForwardPass(rbt,cst,params,xbar,ubar,du,K);
            obj.Update_iter();
            obj.J_pushback(Vbar);
            %%% start iteration
            while obj.iter <= params.Max_iter
                if params.Debug == 1
                    fprintf('[INFO]: Iteration %3d   ||  Cost %.12e \n',obj.iter,Vbar);
                end
                success = 0;
                while success == 0
                    if params.Debug == 1
                        fprintf('\t \t [SubInfo]: Reg=%5f\n',obj.Reg);
                    end
                    [dV,Vx,Vxx,du,K,success] = obj.BackwardPass(rbt,cst,xbar,ubar,params);
                    if success == 0
                        % Need increase Reg factor (min incremental is 1e-3)
                        obj.Reg = max(obj.Reg*2, 1e-3);
%                         obj.Reg = obj.Reg * 2;
                    end   
                end
                Vprev = Vbar;
                % Set regularization back to 0 for next backward pass
                obj.Reg = 0.000; 
                
                %%% Forward Pass
                [Vbar,xbar,ubar] = obj.ForwardIteration(xbar,ubar,Vprev,du,K,dV,rbt,cst,params);
                obj.J_pushback(Vbar);
                obj.Update_iter();
                change = Vprev - Vbar;
                if change < params.stop
                    break
                end
            end
        end
    end
end

