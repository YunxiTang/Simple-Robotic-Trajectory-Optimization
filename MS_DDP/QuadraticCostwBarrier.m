classdef QuadraticCostwBarrier < cost_mdl
    %QUADRATICCOSTWBARRIER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R                   % path input cost (integral{1/2*(u'*R*u*dt)})
        Q                   % path state cost
        Q_f                 % terminal cost (1/2*(x-x_star)'*Qf*(x-x_star))
        u_max               % limit control input
        u_min               % minimum control input
        iter                % iterations
        epi_cost            % costs of all the episodes
    end
    
    methods
        function obj = QuadraticCostwBarrier(Q_f, Q, R, umax,umin)
            %QUADRATICCOSTWBARRIER Construct an instance of this class
            fprintf("[INFO]: Initializing Quadratic Cost Function. \n");
            obj.Q_f = Q_f;
            obj.R = R;
            obj.Q = Q;
            obj.u_max = umax;
            obj.u_min = umin;
            obj.iter = 0;
            obj.epi_cost = [];
        end
        
        function [] = update_iter(obj)
            obj.iter = obj.iter + 1;
        end
        
        function [] = stack_push_J(obj, cost)
            obj.epi_cost = [obj.epi_cost cost];
        end
        
        function [J] = calc_J(obj,x_hst,u_hst,x_ref,u_ref,dt)
            % calculate the cost of single trail
            [S, L, Nu,~] = size(u_hst);
            [~, ~, Nx,~] = size(x_hst);
            J = 0.0;
            for i=1:S
                for j=1:L
                    x = squeeze(x_hst(i, j, :, :));
                    u = squeeze(u_hst(i, j, :, :));
                    x = reshape(x,[Nx,1]);
                    u = reshape(u,[Nu,1]);
                    J = J + obj.L(x, u, x_ref, u_ref, dt);
                end
            end
            obj.stack_push_J(J);
        end
        
        function [] = update_Qf(obj,H)
            obj.Q_f = H;
        end
        
        function phi = phi(obj,x_f,x_ref)
            % phi terminal state cost
            Nx = numel(x_ref);
            x_f = reshape(x_f, Nx, 1);
            phi = 1/2 .* (x_f-x_ref).' * obj.Q_f * (x_f-x_ref);
        end
        
        % terminal state cost derivatives
        function phi_x = phi_x(obj,x_f,x_ref)
            Nx = numel(x_ref);
            x_f = reshape(x_f, Nx, 1);
            phi_x = obj.Q_f * (x_f - x_ref);
        end
        
        function phi_xx = phi_xx(obj, x_f, x_ref)
            phi_xx = obj.Q_f;
        end
        
        % barrier functions
        function B = bar(obj,x,u,dt)
            B = (1 / obj.iter) * ((obj.u_max - u).^(-1) + (u - obj.u_min).^(-1))*dt;
        end
        
        function B_u = bar_u(obj,x,u,dt)
            B_u = (1 / obj.iter) * ((1).*(obj.u_max - u).^(-2) + (-1).*(u - obj.u_min).^(-2))*dt;
        end
        
        function B_uu = bar_uu(obj,x,u,dt)
            B_uu = (1 / obj.iter) * ((2).*(obj.u_max - u).^(-3) + (2).*(u - obj.u_min).^(-3))*dt;
        end
        
        function B_x = bar_x(obj,x,u,dt)
            B_x = zeros(numel(x),1);
        end
        
        function B_xx = bar_xx(obj,x,u,dt)
            B_xx = zeros(numel(x));
        end
        
        % path (or running) cost
        function L = L(obj, x, u, x_ref, u_ref, dt)
            Nx = numel(x_ref);
            x = reshape(x, Nx, 1);
            L = 1/2 .* (u-u_ref).' * obj.R * (u-u_ref) .* dt  + ...
                1/2 .* (x-x_ref).' * obj.Q * (x-x_ref).* dt + obj.bar(x, u, dt);  
        end
        
        % running cost derivatives
        function L_x = L_x(obj, x, u, x_ref, u_ref, dt)
            Nx = numel(x_ref);
            x = reshape(x, Nx, 1);
            L_x = obj.Q * (x - x_ref).* dt + obj.bar_x(x, u, dt);
        end
        
        function L_u = L_u(obj, x, u, x_ref, u_ref, dt)
            Nx = numel(x_ref);
            x = reshape(x, Nx, 1);
            L_u = obj.R * (u-u_ref) .* dt + obj.bar_u(x, u, dt);
        end
        
        function L_xx = L_xx(obj, x, u, x_ref, u_ref, dt)
            Nx = numel(x_ref);
            x = reshape(x, Nx, 1);
            L_xx = obj.Q .* dt + obj.bar_xx(x, u, dt);
        end
        
        function L_uu = L_uu(obj, x, u, x_ref, u_ref, dt)
            L_uu = obj.R .* dt + obj.bar_uu(x, u, dt);
        end
        
        function L_ux = L_ux(obj, x, u, x_ref, u_ref, dt)
            L_ux = zeros(numel(u),numel(x));
        end
        
        function L_xu = L_xu(obj, x, u, x_ref, u_ref, dt)
            L_xu = zeros(numel(x),numel(u));
        end
    end
end

