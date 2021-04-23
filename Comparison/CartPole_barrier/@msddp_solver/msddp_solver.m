classdef msddp_solver < handle
    %MSDDPSOLVER multiple shooting DDP solver
    
    properties
        M,                  % shooting phase
        L,                  % length of shooting phase
        Reg = 1e-3,         % how much to regularize
        Reg_Type = 1,       % 1->reg Quu (Default) / 2->reg Vxx
        eps = 1.0,          % eps: line-search parameter  
        gamma = 0.1,        % threshold to accept a FW step
        beta = 0.5,         % for line-search backtracking
        iter = 0,           % count iterations
        Jstore = [],        % store full costs
        J_real = [],        % store real costs
        Contract_Rate = []  % store constraction rate
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
            % store costs
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
        
        function [] = solver_Callback(obj,xbar,ubar,params)
            % For plot
            nx = params.nx;
            dft = obj.CalDefect(xbar,params);
            for i = 1:obj.M
                clr = [i, 0.5 * i, 0.5 * i] / obj.M;
                figure(111);
                plot(xbar{i}(1,:),xbar{i}(2,:),'Color',clr,'LineWidth',2.0);hold on; 
            end
            axis equal;
            title('Phase Portrait','Interpreter','latex','FontSize',20);
            grid on;
            hold off;
            
            for i = 1:obj.M
                clr = [i, 0.5 * i, 0.5 * i] / obj.M;
                figure(222);
                plot(params.t{i}(1,1:end-1),ubar{i}(:,1:end-1),'Color',clr,'LineWidth',2.0);hold on; 
            end
            title('Control Input','Interpreter','latex','FontSize',20);
            grid on;
            hold off;
            
            figure(555);
            title('Defects','Interpreter','latex','FontSize',15);
            for k=1:nx
                subplot(nx,1,k);
                plot(dft(k,:),'s','LineWidth',2.0);
                grid on;
            end
            
        end
        
        function [J_idx,Jr_idx, xsol,usol] = simulate_phase(obj,rbt,cst,path_constraint,final_constraint,params,idx,x0,xbar,ubar,du,K)
            % simulate each inter-phase
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
                usol(:,i) = ui;
                
                % Sum up costs
                J_idx = J_idx + cst.l_cost(xi, ui, xref, uref) + path_constraint.penalty(xi, ui)*params.dt;
                Jr_idx = Jr_idx + cst.l_cost(xi, ui, xref, uref);
                % Propagate dynamics
                xi = rbt.rk(xi, ui, params.dt);
                if isnan(xi)
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
       
        function [J,Jr,xbar,ubar] = Init_Forward(obj,rbt,cst,path_constraint,final_constraint,params)
            % init forward simulation
            xbar = cell(obj.M, 1);
            ubar = cell(obj.M, 1);
            J = 0;
            Jr = 0;
            if params.warm_start == 1
                x_al = load('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\Comparison\CartPole_ms_constrained\x_al.mat');
                u_al = load('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\Comparison\CartPole_ms_constrained\u_al.mat');
                x_warm = x_al.xbar;
                u_warm = u_al.ubar;
                xbar = x_warm;
                ubar = u_warm;
            else
                for k = 1:obj.M
                    % make initial guess of intermediate states
                    xbar{k} = kron(ones(1, obj.L), params.xf);
                    ubar{k} = zeros(params.nu, obj.L);
                    xbar{k}(:,1) = params.x0 + (k-1) * (params.xf - params.x0) ./ obj.M;
                end
            end
            du = zeros(params.nu, obj.L);
            K  = zeros(params.nu, params.nx, obj.L);
            for i=1:obj.M
                x0 = xbar{i}(:,1);
                [J_idx,Jr_idx,xbar{i},ubar{i}] = obj.simulate_phase(rbt,cst,path_constraint,final_constraint,params,i,x0,xbar{i},ubar{i},du,K);
                J = J + J_idx;
                Jr = Jr + Jr_idx;
            end
            obj.J_pushback(J);
            obj.Jreal_pushback(Jr);
        end
        
        
        function [dft] = CalDefect(obj,xbar,params)
            % compute all the defects
            dft = zeros(params.nx, params.shooting_phase-1);
            for i = 1:obj.M-1
                dft(:,i) = xbar{i}(:,end) - xbar{i+1}(:,1);
            end
        end
        
        function [J,Jr,x,u] = ForwardPass(obj,rbt,cst,path_constraint,final_constraint,params,xbar,ubar,du,K)
            % foward simulation
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
                [J_idx,Jr_idx,x{i},u{i}] = obj.simulate_phase(rbt,cst,path_constraint,final_constraint,params,i,x0,xbar{i},ubar{i},du{i},K{i});
                J = J + J_idx;
                Jr = Jr + Jr_idx;
            end
            
        end
        
        function [dV,Vx,Vxx,du,K,success] = BackwardPass(obj,rbt,cst,path_constraint,final_constraint,xbar,ubar,dft,params)
            % backward propogation
            success = 1;
            xref = params.xf;
            uref = zeros(params.nu, 1);
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
                    [Qx,Qu,Qxx,Quu,Qux,Qxu] = cst.Qcs_info(rbt,cst,x_ij,u_ij,xref,uref,Vx_ij,Vxx_ij,params,path_constraint);
                    % regularization
                    Quu_reg = Quu + eye(params.nu)*obj.Reg;
                    % Make sure Quu is PD, if not, exit and increase regularization
                    [~, FLAG] = chol(Quu_reg-eye(params.nu)*1e-9);
                    if FLAG > 0 
                        % Quu is not PD, then increase Reg factor until Quu
                        % is PD
                        success = 0;
                        if params.Debug == 1
                            fprintf(' \t \t [SubSubInfo]: Reg=%5f:  Regularization <-- Increase Reg. \n', obj.Reg);
                        end
                        return
                    end
                    
                    % Standard Recursive Equations
                    % add boxQP solver here
                    if params.qp == 1 && obj.iter > 10
                        lb = params.umin * ones(params.nu, 1);
                        ub = params.umax * ones(params.nu, 1);
                        lower = lb - u_ij;
                        upper = ub - u_ij;
                        [kff,~,R,free] = boxQP(Quu_reg, Qu, lower, upper);
                        kfb = zeros(params.nu, params.nx);
                        if any(free)
                            Lfree = -R\(R'\Qux(free,:));
                            kfb(free,:) = Lfree;
                        end
                    else
                        % without boxQP 
                        % more numerical stable
                        [R, ~] = chol(Quu_reg);
                        kff = -R\(R'\Qu);
                        kfb = -R\(R'\Qux);
                    end
                    du{i}(:,j) = kff;
                    K{i}(:,:,j) = kfb;
                    
                    Vx{i}(:,j)  = Qx + kfb'*Quu*kff + kfb'*Qu + Qxu*kff ;
                    Vxx{i}(:,:,j) = Qxx + Qxu*kfb + kfb'*Qux + kfb'*Quu*kfb;
                    alpha = obj.eps;
                    delta1 = kff'*Qu + gap'*(Vx{i}(:,j) - Vxx{i}(:,:,j)*x_ij);
                    delta2 = kff'*Quu*kff + gap'*(2*Vxx{i}(:,:,j)*x_ij-Vxx{i}(:,:,j)*gap);
                    dV = alpha * delta1 + 1/2*alpha^2*delta2;
                end
            end
            if params.Debug == 1
                fprintf('\t \t [SubInfo]: Reg=%5f\n',obj.Reg);
            end
        end
        
        function [V,Vr,x,u] = ForwardIteration(obj,xbar,ubar,Vprev,du,K,dV,rbt,cst,path_constraint,final_constraint,params)
            % Check the Forward Pass is accepted or not (IMPORTANT!!!)
            % if not, adjust the line-search factor  (Armijo backtracking
            % line search)
            V = 0;
            Vr = 0;
            obj.eps = 1.0;
            alpha = obj.eps;
            while alpha > 1e-9
                % Try a step
                alpha = obj.eps;
                [V,Vr,x,u] = obj.ForwardPass(rbt,cst,path_constraint,final_constraint,params,xbar,ubar,du,K);
                ratio = (V - Vprev)/(dV);
                if params.Debug == 1
                    fprintf(' \t \t \t Alpha=%.3e \t Actual Reduction=%.3e \t Expected Reduction=%.3e \t Ratio=%.3e\n',...
                          alpha,V-Vprev, obj.gamma*dV,ratio);
                end
                if obj.iter == 0
                    break
                end
                if dV < 0 && V < Vprev + obj.gamma * dV
                    break
                end
                if dV >= 0 && ratio < 5
                    break
                end
                % Else, do backtracking
                obj.eps = obj.beta * obj.eps;
            end
            obj.J_pushback(V);
            obj.Jreal_pushback(Vr);
            obj.Rate_pushback(ratio);
            if params.plot == 1 && mod(obj.iter,2) == 0
               obj.solver_Callback(x,u,params);
            end
        end
        
        function [xsol, usol, Ksol] = Solve(obj,rbt,cst,path_constraint,final_constraint,params)
            % solve OCP
            % init rolling out
            [Vbar,Vr,xbar,ubar] = obj.Init_Forward(rbt,cst,path_constraint,final_constraint,params);
            
            % plot initial trajectories
            if params.plot == 1
                obj.solver_Callback(xbar,ubar,params);
            end
            [dft] = obj.CalDefect(xbar,params);
            while obj.iter <= params.Max_iter
                success = 0;
                % run backward pass
                while success == 0 
                    [dV,Vx,Vxx,du,K,success] = obj.BackwardPass(rbt,cst,path_constraint,final_constraint,xbar,ubar,dft,params);
                    obj.Reg = max(obj.Reg * 2, 1e-3);
                end
                Vprev = Vbar;
                if params.Debug == 1
                    fprintf('[INFO]: Iteration %3d   ||  Cost %.12e \n',obj.iter,Vbar);
                end
                % Set regularization back to 0 for next backward pass
                obj.Reg = 1e-3;
                %%% Forward Pass
                [Vbar,Vr,xbar,ubar] = obj.ForwardIteration(xbar,ubar,Vprev,du,K,dV,rbt,cst,path_constraint,final_constraint,params);
                
                [dft] = obj.CalDefect(xbar,params);
                change = Vprev - Vbar;
                DU = cell2mat(du);
                DU = reshape(DU,(params.nu*params.shooting_phase*obj.L),1);
                if (change) < params.stop && obj.iter > 5 || all(DU<1e-3)
                    break
                end
                obj.Update_iter();
                if mod(obj.iter, 3)==0
                    path_constraint.update_t();
                    final_constraint.update_t();
                end
            end
            [xsol, usol, Ksol] = obj.assemble_solution(xbar, ubar, K, params);
        end
        
        function [xsol, usol, Ksol] = assemble_solution(obj, xbar, ubar, K, params)
            %%% Assemble the final solution
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

