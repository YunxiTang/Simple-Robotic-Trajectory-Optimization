classdef Bipedal < handle
    %BIPEDAL Three Link Bipedal Robot
    
    properties
        Name = 'Bipedal',
        m1 = 7;
        m2 = 7;
        m3 = 15;

        l1 = 0.5;
        l2 = 0.5;
        l3 = 0.35;

        g = 9.81;
    end
    
    methods
        function obj = Bipedal()
            %BIPEDAL 
            fprintf('[INFO]: Creating A Three-link Bipedal Robot.');
        end
        
        function [M, C, G, B] = EoM(obj, x)
            q1 = x(1);
            q2 = x(2);
            q3 = x(3);
            dq1 = x(4);
            dq2 = x(5);
            dq3 = x(6);

            M = eval_M_tmp(obj.l1,obj.l2,obj.l3,obj.m1,obj.m2,obj.m3,q1,q2,q3);
            C = eval_C_tmp(dq1,dq2,dq3,obj.l1,obj.l2,obj.l3,obj.m2,obj.m3,q1,q2,q3);
            G = eval_G_tmp(obj.g,obj.l1,obj.l2,obj.l3,obj.m1,obj.m2,obj.m3,q1,q2,q3);
            B = eval_B_tmp();
        end
        
        function xd = Dynamics(obj, t, x, u)
            %METHOD1 Dynamics
            % From "Feedback Control of Dynamic Bipedal Robot Locomotion"
            [M, C, G, B] = obj.EoM(x);
            dq = [x(4); x(5); x(6)];
            xd = zeros(6, 1);
            xd(1) = x(4);
            xd(2) = x(5);
            xd(3) = x(6);
            xd(4:6) = M \ (-C*dq - G + B*u);
        end
        
        function x_next = rk(obj, x, u, dt)
           k1 = obj.Dynamics(0, x,             u);
           k2 = obj.Dynamics(0, x + 0.5*dt*k1, u);
           k3 = obj.Dynamics(0, x + 0.5*dt*k2, u);
           k4 = obj.Dynamics(0, x +     dt*k3, u);
           x_next = x + dt/6*(k1+2*k2+2*k3+k4);
        end
        
        function [fx,fu] = getLinSys(obj,qbar,ubar,dt)
            dh = 1e-4;
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
    methods (Static)
        M = eval_M_tmp(l1,l2,l3,m1,m2,m3,q1,q2,q3);
        C = eval_C_tmp(dq1,dq2,dq3,l1,l2,l3,m2,m3,q1,q2,q3);
        G = eval_G_tmp(g,l1,l2,l3,m1,m2,m3,q1,q2,q3);
        B = eval_B_tmp();
        T = eval_T_tmp(dq1,dq2,dq3,l1,l2,l3,m1,m2,m3,q1,q2,q3);
        V = eval_V_tmp(dq1,dq2,dq3,g,l1,l2,l3,m1,m2,m3,q1,q2,q3);
    end
end

