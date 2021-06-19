classdef planararm < handle
    
    properties
        Name = 'Planar Arm',
        NJ,
        rbt,
        Nx,
        Nu,
        B
    end
    
    methods
        function obj = planararm(mdl, nj)
            % CartPole Construct an instance of this class
            obj.rbt = mdl;
            obj.NJ  = nj;
            obj.Nx = 2 * nj;
            obj.Nu = nj;
            obj.B  = diag([1; zeros(nj-1, 1)]);
            disp('[INFO]: Creating A planar arm.');
        end
        
        function xd = Dynamics(obj, t, x, u)
            %%% planar arm dynamics (ABA FD)
            nj = obj.NJ;
            q = x(1:nj, :);
            qd = x((nj+1):2*nj, :);
            qdd = FDab(obj.rbt, q, qd, u);
            xd = [qd;qdd];
        end
        
        function q_next = rk(obj, q, u, dt)
           k1 = obj.Dynamics(0, q,             u);
           k2 = obj.Dynamics(0, q + 0.5*dt*k1, u);
           k3 = obj.Dynamics(0, q + 0.5*dt*k2, u);
           k4 = obj.Dynamics(0, q +     dt*k3, u);
           q_next = q + dt/6*(k1+2*k2+2*k3+k4);
        end
        
        function [fx,fu] = getLinSys(obj,qbar,ubar,dt)
            dh = 1e-6;
            Jx = zeros(obj.Nx,obj.Nx);
            Ju = zeros(obj.Nx,obj.Nu);
            Hx = eye(obj.Nx) * dh;
            Hu = eye(obj.Nu) * dh;
            for i=1:obj.Nx
                Jx(:,i) = (obj.rk(qbar + Hx(:,i), ubar, dt) - obj.rk(qbar - Hx(:,i), ubar, dt)) / (2*dh);
            end

            for i=1:obj.Nu
                Ju(:,i) = (obj.rk(qbar, ubar+ Hu(:,i), dt) - obj.rk(qbar , ubar - Hu(:,i), dt)) / (2*dh);
            end
            fx = Jx;
            fu = Ju;
        end
        
    end
%     methods (Static)
%         [fx,fu] = getLinSys(in1,in2,internal1);
%     end
end

