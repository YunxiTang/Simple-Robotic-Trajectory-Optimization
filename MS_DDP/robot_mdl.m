classdef (Abstract) robot_mdl < handle
    %ROBOT_MODEL an abstract class defined for robot
    
    properties
        
    end
    
    methods
        % EoM
        [M,C,G,F,B] = Sim_EoM(obj,x)
        
        % dynamical equation
        dxdt = Dynamics(obj, ~, x, u);
        
        % linearized system equations with dt involved
        [A, B] = getLinSys(obj, x, u);
        [f_x, f_u] = f_d(obj, x, u, dt);
        
        % forward kinematics
        x_next = rk45(obj, x, u, dt);
        
        % run animation
        [] = run_animation(obj, t_data, state_data);
    end
end
