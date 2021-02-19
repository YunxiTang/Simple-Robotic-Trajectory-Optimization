classdef rbt_mdl
    %RBT_MODEL robot model
    
    properties
        L = 1.0,
        m = 1.0,
        r = 0.5;
        I = 0.33;
        b = 0.12;
        g = 9.81
    end
    
    methods
        function obj = rbt_mdl()
            %RBT_MODEL Construct an instance of this class
            disp('[INFO]: Creating Pendulum Model.');
        end
        
        function [M, C, G, F, B] = EoM(obj, x)
            q  = x(1,:);
            M = [obj.m*obj.r^2+obj.I];
            C = [0];
            G = [obj.m*obj.g*obj.r*sin(q)];
            F = [obj.b];
            B = [1];
        end
        
        function xd = Dynamics(obj, t, x, u)
            qd = x(2,:);
            [M,C,G,F,B] = obj.EoM(x);
            xd = [qd; 
                  M \ (B*u-C*qd-F*qd-G)];
        end
    end
    
    methods (Static)
        [fx,fu] = getLinSys(in1,u1);
    end
end

