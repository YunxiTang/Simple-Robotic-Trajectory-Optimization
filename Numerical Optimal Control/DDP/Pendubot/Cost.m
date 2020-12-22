classdef (Abstract) Cost < handle
    %COST An abstract cost class for terminal and path cost 
    
    methods
        % terminal state cost
        phi = phi(obj, x_f, x_star);
        
        % terminal state cost derivatives
        phi_x = phi_x(obj, x_f, x_star);
        phi_xx = phi_xx(obj, x_f, x_star);
        
        % running cost
        L = L(obj, x, u, dt);
        
        % running cost derivatives
        L_x = L_x(obj, x, u, dt);
        L_u = L_u(obj, x, u, dt);
        L_xx = L_xx(obj, x, u, dt);
        L_uu = L_uu(obj, x, u, dt);
        L_ux = L_ux(obj, x, u, dt);
        L_xu = L_xu(obj, x, u, dt);
    end
end

