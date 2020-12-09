classdef QuadraticCost < Cost
    %QUADRATICCOST Quadratic terminal and path cost which inherits from
    %abstract Cost class (Qf and R are needed)
    
    properties
        R                   % path cost (integral{1/2*(x'*R*x*dt)})
        Q_f                 % terminal cost (1/2*(x-x_star)'*Qf*(x-x_star))
    end
    
    methods
        function obj = QuadraticCost(Q_f,R)
            %QUADRATICCOST Construct an instance of this class
            %   Detailed explanation goes here
            obj.Q_f = Q_f;
            obj.R = R;
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
        
        % path (or running) cost
        function L = L(obj, x, u, dt)
            L = 1/2 .* u.' * obj.R * u .* dt;  
        end
        
        % running cost derivatives
        function L_x = L_x(obj, x, u, dt)
            L_x = zeros(numel(x),1);
        end
        
        function L_u = L_u(obj, x, u, dt)
            L_u = obj.R * u .* dt;
        end
        
        function L_xx = L_xx(obj, x, u, dt)
            L_xx = zeros(numel(x));
        end
        
        function L_uu = L_uu(obj, x, u, dt)
            L_uu = obj.R .* dt;
        end
        
        function L_ux = L_ux(obj, x, u, dt)
            L_ux = zeros(numel(u),numel(x));
        end
        
        function L_xu = L_xu(obj, x, u, dt)
            L_xu = zeros(numel(x),numel(u));
        end
            
    end
end

