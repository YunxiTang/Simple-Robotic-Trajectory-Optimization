classdef Car < handle
    %Car A 2D Car model
    
    properties
        Name = '2D Car'
        Nx = 4,
        Nu = 2
    end
    
    methods
        function obj = Car()
            %Car Construct an instance of this class
            disp('[INFO]: Creating 2D Car Model.');
        end
        
        function [dq] = Dynamics(obj,t,q,u)
                %DYNAMICS 
                % dynamics

                x = q(1);
                y = q(2);
                theta = q(3);
                v = q(4);

                u_theta = u(1);
                u_v = u(2);

                dq = [v*sin(theta);
                      v*cos(theta);
                      u_theta * v;
                      u_v]; 
        end
        
        function q_next = rk(obj, q, u, dt)
            x = q(1);
            y = q(2);
            theta = q(3);
            v = q(4);
            u_the = u(1);
            u_v = u(2);
            x_next = x + dt * v * sin(theta);
            y_next = y + dt * v * cos(theta);
            theta_next = theta + dt * u_the * v;
            v_next = v + dt * u_v;
            q_next = [x_next; y_next; theta_next; v_next];
        end
    end
    methods (Static)
        [fx,fu] = getLinSys(in1,in2,dt);
    end
end

