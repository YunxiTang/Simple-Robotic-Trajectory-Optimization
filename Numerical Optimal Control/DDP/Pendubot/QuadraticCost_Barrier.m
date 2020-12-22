classdef QuadraticCost_Barrier < Cost
    %QUADRATICCOST Quadratic terminal and path cost which inherits from
    %abstract Cost class (Qf and R are needed)
    
    properties
        R                   % path cost (integral{1/2*(u'*R*u*dt)})
        Q
        Q_f                 % terminal cost (1/2*(x-x_star)'*Qf*(x-x_star))
        u_max               % limit control input
        u_min               % minimum control input
        iter                % iteration time
    end
    
    methods
        function obj = QuadraticCost_Barrier(Q_f,Q,R, umax,umin)
            %QUADRATICCOST Construct an instance of this class
            %   Detailed explanation goes here
            obj.Q_f = Q_f;
            obj.R = R;
            obj.Q = Q;
            obj.u_max = umax;
            obj.u_min = umin;
            obj.iter = 0;
        end
        
        function [] = update_iter(obj)
            obj.iter = obj.iter + 1;
        end
        
        function phi = phi(obj,x_f,x_star)
            %phi terminal state cost
            phi = 1/2 .* (x_f-x_star).' * obj.Q_f * (x_f-x_star);
        end
        
        % terminal state cost derivatives
        function phi_x = phi_x(obj,x_f,x_star)
            phi_x = obj.Q_f * (x_f - x_star);
        end
        
        function phi_xx = phi_xx(obj, x_f, x_star)
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
        function L = L(obj, x, u, x_star, dt)
            L = 1/2 .* u.' * obj.R * u .* dt + obj.bar(x, u, dt) + ...
                1/2 .* (x-x_star).' * obj.Q * (x-x_star).* dt;  
        end
        
        % running cost derivatives
        function L_x = L_x(obj, x, u, x_star, dt)
            L_x = obj.Q * (x - x_star).* dt + obj.bar_x(x, u, dt);
        end
        
        function L_u = L_u(obj, x, u, x_star, dt)
            L_u = obj.R * u .* dt + obj.bar_u(x, u, dt);
        end
        
        function L_xx = L_xx(obj, x, u, x_star, dt)
            L_xx = obj.Q.* dt + obj.bar_xx(x, u, dt);
        end
        
        function L_uu = L_uu(obj, x, u, x_star, dt)
            L_uu = obj.R .* dt + obj.bar_uu(x, u, dt);
        end
        
        function L_ux = L_ux(obj, x, u, x_star, dt)
            L_ux = zeros(numel(u),numel(x));
        end
        
        function L_xu = L_xu(obj, x, u, x_star, dt)
            L_xu = zeros(numel(x),numel(u));
        end
            
    end
end

