classdef constraint < handle
    %CONSTRAINT add constraint here, solve with Augmented Lagrange Method
    %%% State-Control Constraint
        % cineq(x, u) <= 0
    
    properties
        n_ineq = 0.0,       % Number of inequality constraint
        x_obstacle = [],    % [nx+1, n_ineq], store obstacles
        active_set = [],     % check the active constraints
        dt
    end
    
    methods
        function obj = constraint(obstacles, dt)
            %CONSTRAINT Construct an instance of constraint
            disp("[INFO]: Creating Constraint Model.");
            [~, num_ineq] = size(obstacles);
            obj.n_ineq = num_ineq;
            obj.x_obstacle = obstacles;
            obj.dt = dt; 
        end
        
        function Imu = active_check(obj, x, u, lambda, mu)
            %METHOD1 check the active set
            Imu = zeros(obj.n_ineq, obj.n_ineq);
            h = obj.c_ineq(x, u);
            for i=1:obj.n_ineq
                if h(i) > -1e-5 || lambda(i) > 0.0
                    Imu(i,i) = mu(i);
                end
            end
        end
        
        function h = c_ineq(obj, x, u)
            % ineqaulity constraint (make sure it is less/equal (h<=0) than 0)
            % pick out the constrained subspace
            H = [1 0 0 0 0 0;
                 0 1 0 0 0 0];
            
            for k=1:obj.n_ineq
                r = obj.x_obstacle(3,k);
                x_obs = obj.x_obstacle(1:2,k);
                h(k) = r*r - (H*x - x_obs)' * (H*x - x_obs);
            end
            h = reshape(h,[obj.n_ineq,1]);
        end
        
        function [AL_J] = AL_Term(obj, x, u, lambda, Imu)
            AL_j = (lambda + 1/2 * Imu * obj.c_ineq(x,u))' * obj.c_ineq(x,u);
            AL_J = AL_j * obj.dt;
        end
      
        
        function Flag = plot_obstacle(obj, map)
            if nargin > 1
                figure(map);
            end
            aplha=0:0.01:2*pi;
            for i=1:obj.n_ineq
                % the i-th circle
                r = obj.x_obstacle(3,i); 
                xi_c = obj.x_obstacle(1,i);
                yi_c = obj.x_obstacle(2,i);
                cx = xi_c + r*cos(aplha);
                cy = yi_c + r*sin(aplha);
                clr = abs([sin(i/1.2) sin(i/1.5) sin(i)]);
                plot(cx, cy, '-','Color',clr); hold on;
                fill(cx, cy, clr);
                axis equal;
            end
            xlabel('$x$','Interpreter','latex','FontSize',15);
            ylabel('$y$','Interpreter','latex','FontSize',15);
            title("$2D\;Map$", 'Interpreter','latex','FontSize',15);
            grid on;
            Flag = 1;
        end
    end
    methods (Static)
       %%% After code generation to class, Add the function declaration
       %%% here !!!
       [cx,cu] = algrad(x,u);
    end
end


