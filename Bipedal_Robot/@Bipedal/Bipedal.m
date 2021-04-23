classdef Bipedal < handle
    %BIPEDAL Three Link Bipedal Robot
    
    properties
        Name = 'Bipedal',
        l = 0.5,      % torso length
        r = 1.0,      % leg length
        M_T = 10,    % torso mass
        M_H = 15,    % hip mass
        m = 5,      % leg mass
        g = 9.81      % gravity
    end
    
    methods
        function obj = Bipedal()
            %BIPEDAL 
            fprintf('[INFO]: Creating A Bipedal Robot.');
        end
        
        function qd = Dynamics(obj, t, q, u)
            %METHOD1 Dynamics
            % From "Feedback Control of Dynamic Bipedal Robot Locomotion"
            theta_1 = q(1);
            theta_2 = q(2);
            theta_3 = q(3);
            dtheta_1 = q(4);
            dtheta_2 = q(5);
            dtheta_3 = q(6);
            
            M(1,1) = (5 * obj.m / 4 + obj.M_H + obj.M_T) * obj.r * obj.r;
            M(1,2) = -0.5 * obj.m * obj.r * obj.r * cos(theta_1 - theta_2);
            M(1,3) = obj.M_T * obj.r * obj.l * cos(theta_1 - theta_3);
            M(2,1) = M(1,2);
            M(2,2) = 0.25 * obj.m * obj.r * obj.r;
            M(2,3) = 0;
            M(3,1) = M(1,3);
            M(3,2) = M(2,3);
            M(3,3) = obj.M_T * obj.l * obj.l;
            
            C(1,1) = 0.0;
            C(1,2) = -0.5 * obj.m * obj.r * obj.r * sin(theta_1 - theta_2) * dtheta_2;
            
        end
    end
end

