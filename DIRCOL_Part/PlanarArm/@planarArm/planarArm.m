classdef planarArm
    %PLANARARM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name = 'Planar Arm'
        Nj,
        mdl,
        Nx,
        Nu,
        B
    end
    
    methods
        function obj = planarArm(rbt, nj)
            %PLANARARM Construct an instance of this class
            %   import robot model
            obj.mdl = rbt;
            obj.Nj = nj;
            obj.Nx = 2 * nj;
            obj.Nu = nj;
            obj.B  = diag([1; zeros(nj-1, 1)]);
        end
        
        function xd = Dynamics(obj, t, x, u)
            %%% planar arm dynamics (ABA FD)
            nj = obj.Nj;
            q = x(1:nj, :);
            qd = x((nj+1):2*nj, :);
            qdd = FDab(obj.mdl, q, qd, obj.B *u);
            xd = [qd;qdd];
        end
    end
end

