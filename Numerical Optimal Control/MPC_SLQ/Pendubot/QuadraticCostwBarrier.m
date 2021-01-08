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
    end
    
    methods
        function obj = QuadraticCostwBarrier(Q_f, Q, R, umax,umin)
            %QUADRATICCOSTWBARRIER Construct an instance of this class
            obj.Q_f = 500*Q_f;
            obj.R = R;
            obj.Q = Q;
            obj.u_max = umax;
            obj.u_min = umin;
            obj.iter = 1;
        end
        
        function [] = update_iter(obj)
            obj.iter = obj.iter + 1;
        end
        
        function [] = update_Qf(obj,H)
            obj.Q_f = H;
        end
        
        function phi = phi(obj,x_f,x_ref)
            % phi terminal state cost
            phi = 1/2 .* (x_f-x_ref).' * obj.Q_f * (x_f-x_ref);
        end
        
        % terminal state cost derivatives
        function phi_x = phi_x(obj,x_f,x_ref)
            phi_x = obj.Q_f * (x_f - x_ref);
        end
        
        function phi_xx = phi_xx(obj, x_f, x_ref)
            phi_xx = obj.Q_f;
        end
        
        % path (or running) cost
        function L = L(obj, x, u, x_ref, u_ref, dt)
            L = 1/2 .* (u-u_ref).' * obj.R * (u-u_ref) .* dt  + ...
                1/2 .* (x-x_ref).' * obj.Q * (x-x_ref).* dt;  
        end
        
        % running cost derivatives
        function L_x = L_x(obj, x, u, x_ref, u_ref, dt)
            L_x = obj.Q * (x - x_ref).* dt ;
        end
        
        function L_u = L_u(obj, x, u, x_ref, u_ref, dt)
            L_u = obj.R * (u-u_ref) .* dt ;
        end
        
        function L_xx = L_xx(obj, x, u, x_ref, u_ref, dt)
            L_xx = obj.Q .* dt ;
        end
        
        function L_uu = L_uu(obj, x, u, x_ref, u_ref, dt)
            L_uu = obj.R .* dt ;
        end
        
        function L_ux = L_ux(obj, x, u, x_ref, u_ref, dt)
            L_ux = zeros(numel(u),numel(x));
        end
        
        function L_xu = L_xu(obj, x, u, x_ref, u_ref, dt)
            L_xu = zeros(numel(x),numel(u));
        end
    end
end

