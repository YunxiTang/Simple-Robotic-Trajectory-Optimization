classdef RR_robot < handle
    %Car A 2D Car model
    
    properties
        Name = 'RR_Robot',
        l1 = 1.0,
        l2 = 1.0,
        m1 = 1.0,
        m2 = 1.0,
        lc1 = 0.5,
        lc2 = 0.5,
        
        b1 = 0.01,
        b2 = 0.01,
        
        I1 = 0.33,
        I2 = 0.33,
        
        g = 9.81
    end
    
    methods
        function obj = RR_robot()
            %Car Construct an instance of this class
            disp('[INFO]: Creating RR_Robot Model.');
        end
        
        function [M,C,G,F,B] = EoM(obj,x)
            % INPUTS:
            %    model: struct
            %    x: [4,1] = [q1; q2; q1d; q2d]
            %
            % OUTPUTS:
            %    M: [2,2] = inertia matrix
            %    C: [2,2] = coriolis and centrifugal term
            %    G: [2,1] = gravitational term
            %    F: [2,1] = Fiction force term
                q1 = x(1); q2 = x(2);
                dq1 = x(3); dq2 = x(4);

                M11 = obj.I1 + obj.I2 + obj.m2*obj.l1*obj.l1 + 2*obj.m2*obj.l1*obj.lc2*cos(q2);
                M12 = obj.I2 + obj.m2*obj.l1*obj.lc2*cos(q2);
                M21 = M12; M22 = obj.I2;
                M = 1.0*[M11 M12;
                     M21 M22];

                C11 = -2*obj.m2*obj.l1*obj.lc2*sin(q2)*dq2;
                C12 = -obj.m2*obj.l1*obj.lc2*sin(q2)*dq2;
                C21 = obj.m2*obj.l1*obj.lc2*sin(q2)*dq1;
                C22 = 0;
                C = 1.0*[C11 C12;
                         C21 C22];

                G1 = obj.m1*obj.g*obj.lc1*sin(q1)+obj.m1*obj.g*(obj.l1*sin(q1)+obj.lc2*sin(q1+q2));
                G2 = obj.m2*obj.g*obj.lc2*sin(q1+q2);
                G = [G1;G2];
                F = 1.0*[obj.b1 0;
                         0 obj.b2];
                B = [1 0;0 1];
        end
        
        function dxdt = Dynamics(obj,t,x,u)
            % nonlinear dynamics
            dq = x(3:4);
            [M,C,G,F,B] = obj.EoM(x);

            ddq = M \ (B*u - G - C*dq - F*dq); 
            dxdt = [dq;
                    ddq];
        end
        
        function x_next = rk(obj, x, u, dt)
            k1 = obj.Dynamics(0, x,             u);
            k2 = obj.Dynamics(0, x + 0.5*dt*k1, u);
            k3 = obj.Dynamics(0, x + 0.5*dt*k2, u);
            k4 = obj.Dynamics(0, x +     dt*k3, u);
            x_next = x + dt/6*(k1+2*k2+2*k3+k4);
        end
        
        function [fx,fu] = getLinSys(obj,qbar,ubar,dt)
            dh = 1e-3;
            Nx = numel(qbar);
            Nu = numel(ubar);
            Jx = zeros(Nx,Nx);
            Ju = zeros(Nx,Nu);
            Hx = eye(Nx) * dh;
            Hu = eye(Nu) * dh;
            for i=1:Nx
                Jx(:,i) = (obj.rk(qbar + Hx(:,i), ubar, dt) - obj.rk(qbar - Hx(:,i), ubar, dt)) / (2*dh);
            end

            for i=1:Nu
                Ju(:,i) = (obj.rk(qbar, ubar+ Hu(:,i), dt) - obj.rk(qbar , ubar - Hu(:,i), dt)) / (2*dh);
            end
            fx = Jx;
            fu = Ju;
        end
    end
%     methods (Static)
%         [fx,fu] = getLinSys(in1,in2);
%     end
end

