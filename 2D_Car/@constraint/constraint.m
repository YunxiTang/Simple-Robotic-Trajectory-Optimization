classdef constraint < handle
    %CONSTRAINT add constraint here, solve with Augmented Lagrange Method
    %%% State Constraint
        % ceq(x) = 0
        % cin(x) <= 0
    %%% Input Constraint
        % umin <= u <= umax
    
    properties
        n_eq,         % number of equality constraint
        n_in,         % number of inequality constraint
        lambda,       % Lagrangian Multiplier for equality constraints
        mu,           % Lagrangian Multiplier for inequality constraints
        rho = 10,     % Multiplier for Augmented Penalty
        phi = 10      % scaling parameter for rho, rho <- phi * rho
    end
    
    methods
        function obj = constraint()
            %CONSTRAINT Construct an instance of constraint
            disp("[INFO]: Creating Constraint Model.");
            obj.n_eq = 
            obj.lambda =;
        end
        
        function outputArg = active_set(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

