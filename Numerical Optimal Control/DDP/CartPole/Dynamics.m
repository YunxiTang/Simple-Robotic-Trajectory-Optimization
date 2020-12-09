classdef (Abstract) Dynamics < handle
    %DYNAMICS An abstract class defined for system dynamics in form of
    %dxdt=f(x,u)
    
    methods
        % Equations of motion in state space representation
        dxdt = F(obj, x, u);
        
        % linearized equations of motion with dt involved
        Phi = Phi(obj,x,u,dt);
        beta = beta(obj,x,u,dt);
    end
end

