classdef final_constraint < handle
    %FINAL_CONSTRAINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_ineq = 0,         % Number of inequality constraint
        final_cons_func,    % path constraints function handle
        active_set = [],    % check the active constraints
    end
    
    methods
        function obj = final_constraint(final_cons, N_ineq)
            %CONSTRAINT Construct an instance of constraint
            disp("[INFO]: Creating Final Constraint Model.");
            obj.final_cons_func = final_cons;
            obj.n_ineq = N_ineq;
        end
        
        function h = c(obj, x)
            % ineqaulity constraint (make sure it is less/equal (h<=0) than 0)
            h = obj.final_cons_func(x);
            h = reshape(h,[obj.n_ineq,1]);
        end
        
        function Imu = active_check(obj, x, lambda, mu)
            %check the active set
            Imu = zeros(obj.n_ineq, obj.n_ineq);
            h = obj.c(x);
            for i=1:obj.n_ineq
                if h(i) > -1e-5 || lambda(i) > 0.0
                    Imu(i,i) = mu(i);
                end
            end
        end
        
        function [AL_J] = AL_Term(obj, x, lambda, Imu)
            AL_j = (lambda + 1/2 * Imu * obj.c(x))' * obj.c(x);
            AL_J = AL_j;
        end

        function M = constraint_violation(obj,x,map)
            if nargin > 2
                figure(map);
            end
            h = obj.c(x);
            M = max(h,[],'all');
        end
        
        function [cx] = algrad(obj,x)
            cx = obj.cx(x);
        end
    end
    methods (Static)
       %%% After code generation to class, Add the function declaration
       %%% here !!!
       cx = cx(in1,u1);
    end
end

