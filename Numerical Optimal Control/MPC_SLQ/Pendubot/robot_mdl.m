classdef (Abstract) robot_mdl < handle
    %ROBOT_MODEL an abstract class defined for robot
    
    properties
        
    end
    
    methods
        % EoM
        [M,C,G,F,B] = EoM(obj,x)
        
        % dynamical equation
        dxdt = f(obj, ~, x, u);
        
        % linearized system equations with dt involved
        [A, B] = getLinSys(obj, x, u);
        [f_x, f_u] = f_x(obj, x, u, dt);
        
        % forward kinematics
        [] = fk(obj, state_data);
        
        % run animation
        [] = run_animation(obj, t_data, state_data);
    end
end

