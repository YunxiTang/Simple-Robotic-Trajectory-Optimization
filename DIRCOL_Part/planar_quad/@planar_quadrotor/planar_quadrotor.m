classdef planar_quadrotor < handle
    %Car A 2D Car model
    
    properties
        Name = 'Planar Quadrotor',
        m = 0.8,
        g = 9.81,
        l = 0.3,
        J = 1.0,
        Nx = 6,
        Nu = 2
    end
    
    methods
        function obj = planar_quadrotor()
            %Car Construct an instance of this class
            disp('[INFO]: Creating Planar Quadrotor Model.');
            obj.J = 0.2*obj.m*obj.l*obj.l;
        end
        
        function xd = Dynamics(obj, t, q, u)
            %%% planar quadrotor dynamics
            theta = q(3);
            
            ddx = -(1 / obj.m) * (u(1) + u(2)) * sin(theta);
            ddy =  (1 / obj.m) * (u(1) + u(2)) * cos(theta) - obj.g;
            ddtheta =  (1 / obj.J) * (obj.l / 2) * (u(2) - u(1));
            xd = [q(4:6);ddx;ddy;ddtheta];
        end
        
        function q_next = rk(obj, q, u, dt)
           k1 = obj.Dynamics(0, q,             u);
           k2 = obj.Dynamics(0, q + 0.5*dt*k1, u);
           k3 = obj.Dynamics(0, q + 0.5*dt*k2, u);
           k4 = obj.Dynamics(0, q +     dt*k3, u);
           q_next = q + dt/6*(k1+2*k2+2*k3+k4);
        end
        
        function [] = animation(obj, T, X, k, Num_Fig)
            N = numel(T);
            fig = figure(Num_Fig);
            
            for i=1:2*k:N
%                 clf(fig);
                xi = X(1,i);
                yi = X(2,i);
                thetai = X(3,i);
                obj.draw_rotor(xi,yi,thetai);
                xlim([0,10]);
                axis equal;
                grid on;
                pause(0.05);
            end
        end
        
        function [] = draw_rotor(obj, x, y, theta)
            xc = x;
            yc = y;
            d = obj.l / 2;
            rotor_L_x = xc - d * cos(theta);
            rotor_R_x = xc + d * cos(theta);
            rotor_L_y = yc - d * sin(theta);
            rotor_R_y = yc + d * sin(theta);
            line([rotor_L_x, rotor_R_x], [rotor_L_y, rotor_R_y],'Color','k','LineWidth',2.0); hold on;
            scatter(xc, yc,'LineWidth',1.5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
            scatter(rotor_L_x, rotor_L_y,'LineWidth',2.5); hold on;
            scatter(rotor_R_x, rotor_R_y,'LineWidth',2.5); hold on;
            
        end
    end
    methods (Static)
        [fx,fu] = getLinSys(in1,in2);
    end
end

