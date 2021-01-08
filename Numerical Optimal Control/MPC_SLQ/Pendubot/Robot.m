classdef (Abstract) Robot < handle
    %DYNAMICS An abstract class defined for system dynamics in form of
    %dxdt=f(x,u)
    
    methods
        % Equations of motion in state space representation
        dxdt = Dynamics(obj, t, x, u);
        
        % linearized equations of motion with dt involved
        [A,B] = getLinSys(obj,x,u);
        Phi = Phi(obj,x,u,dt);
        beta = beta(obj,x,u,dt);
    end
end

