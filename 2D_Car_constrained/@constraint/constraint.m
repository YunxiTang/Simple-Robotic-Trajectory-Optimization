classdef constraint < handle
    %CONSTRAINT add constraint here, solve with Augmented Lagrange Method
    %%% State-Control Constraint
        % cineq(x, u) <= 0
    
    properties
        n_ineq = 0,         % number of inequality constraint
        lambda = 0.0,       % Lagrangian Multiplier for inequality constraints
        rho = 10,           % Multiplier for Augmented Penalty
        phi = 10            % scaling parameter for rho, rho <- phi * rho
        x_obstacle = []     % [nx+1, n_ineq], store obstacles
    end
    
    methods
        function obj = constraint(obstacles)
            %CONSTRAINT Construct an instance of constraint
            disp("[INFO]: Creating Constraint Model.");
            [~, num_ineq] = size(obstacles);
            obj.n_ineq = num_ineq;
            obj.x_obstacle = obstacles;
        end
        
        function outputArg = active_set(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        function g = c_eq(obj, x, u)
            % eqaulity constraint
            g = zeros(obj.n_eq, 1);
            
        end
        
        function h = c_ineq(obj, x, u, x_obs1, r1, x_obs2, r2)
            % ineqaulity constraint (make sure it is less/equal than 0)
            h = zeros(obj.n_ineq, 1);
            H = diag([1.0 1.0 0.0 0.0]);
            h(1) = r1*r1 - (H*x - x_obs1)' * (H*x - x_obs1);
            h(2) = r2*r2 - (H*x - x_obs2)' * (H*x - x_obs2);
        end
        
        function Flag = plot_obstacle(obj, map)
            figure(map);
            aplha=0:pi/40:2*pi;
            for i=1:obj.n_ineq
                % the i-th circle
                r = obj.x_obstacle(3,i); 
                xi_c = obj.x_obstacle(1,i);
                yi_c = obj.x_obstacle(2,i);
                cx = xi_c + r*cos(aplha);
                cy = yi_c + r*sin(aplha);
                clr = abs([sin(i/1.2) sin(i/1.5) sin(i)]);
                plot(cx, cy, '-'); hold on;
                fill(cx, cy, clr);
                axis equal;
            end
            xlabel('$x$','Interpreter','latex','FontSize',15);
            ylabel('$y$','Interpreter','latex','FontSize',15);
            grid on;
            Flag = 1;
        end
    end
end

