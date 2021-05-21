classdef rbt_mdl < handle
    %Car A 2D Car model
    
    properties
        Name = 'CartPole',
        m1 = 2.0, % mass of cart
        m2 = 0.5, % mass of pole
        l = 0.5,  % length of pole
        g = 9.81
    end
    
    methods
        function obj = rbt_mdl()
            % CartPole Construct an instance of this class
            disp('[INFO]: Creating CartPole Model.');
        end
        
        function qd = Dynamics(obj, t, q, u)
            x = q(1);
            theta = q(2);
            dx = q(3);
            dtheta = q(4);
            
            M = [cos(theta)         obj.l;
                 obj.m1 + obj.m2    obj.m2*obj.l*cos(theta)];
            G = [-obj.g * sin(theta);
                  u + obj.m2 * obj.l * dtheta.^2 * sin(theta)];
            acc = M \ G;
            qd = [dx;
                  dtheta;
                  acc]; 
        end
        
        function q_next = rk(obj, q, u, dt)
           k1 = obj.Dynamics(0, q,             u);
           k2 = obj.Dynamics(0, q + 0.5*dt*k1, u);
           k3 = obj.Dynamics(0, q + 0.5*dt*k2, u);
           k4 = obj.Dynamics(0, q +     dt*k3, u);
           q_next = q + dt/6*(k1+2*k2+2*k3+k4);
        end
        
    end
    methods (Static)
        [fx,fu] = getLinSys(in1,in2,dt);
    end
end

