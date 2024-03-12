classdef path_constraint < handle
    %CONSTRAINT add constraint here, solve with Augmented Lagrange Method
    %%% State-Control Constraint, path constraint
        % cineq(x, u) <= 0
    
    properties
        n_ineq = 0,         % Number of inequality constraint
        path_cons_func,     % path constraints function handle
        active_set = [],    % check the active constraints
        dt
    end
    
    methods
        function obj = path_constraint(path_cons, N_ineq, dt)
            %CONSTRAINT Construct an instance of constraint
            disp("[INFO]: Creating Path Constraint Model.");
            obj.path_cons_func = path_cons;
            obj.n_ineq = N_ineq;
            obj.dt = dt; 
        end
        
        function h = c(obj, x, u)
            % ineqaulity constraint (make sure it is less/equal (h<=0) than 0)
            h = obj.path_cons_func(x, u);
            h = reshape(h,[obj.n_ineq,1]);
        end
        
        function Imu = active_check(obj, x, u, lambda, mu)
            %check the active set
            Imu = zeros(obj.n_ineq, obj.n_ineq);
            h = obj.c(x, u);
            for i=1:obj.n_ineq
                if h(i) > -1e-5 || lambda(i) > 0.0
                    Imu(i,i) = mu(i);
                end
            end
        end
        
        function [AL_J] = AL_Term(obj, x, u, lambda, Imu)
            AL_j = (lambda + 1/2 * Imu * obj.c(x,u))' * obj.c(x,u);
            AL_J = AL_j * obj.dt;
        end

        function M = constraint_violation(obj,x,u,map)
            if nargin > 2
                figure(map);
            end
            h = obj.c(x, u);
            M = max(h,[],'all');
        end
        
        function [cx,cu] = algrad(obj,x,u)
            cu = obj.cu(x,u);
            cx = obj.cx(x,u);
        end
    end
    methods (Static)
       %%% After code generation to class, Add the function declaration
       %%% here !!!
       cu = cu(in1,u1);
       cx = cx(in1,u1);
    end
end


