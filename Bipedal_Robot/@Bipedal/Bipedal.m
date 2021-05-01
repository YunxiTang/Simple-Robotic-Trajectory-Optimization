classdef Bipedal < handle
    %BIPEDAL Three Link Bipedal Robot
    
    properties
        Name = 'Bipedal',
        l = 0.5,      % torso length
        r = 1.0,      % leg length
        M_T = 10,     % torso mass
        M_H = 15,     % hip mass
        m = 5,        % leg mass
        g = 9.81      % gravity
    end
    
    methods
        function obj = Bipedal()
            %BIPEDAL 
            fprintf('[INFO]: Creating A Three-link Bipedal Robot.');
        end
        
        function xd = Dynamics(obj, t, x, u)
            %METHOD1 Dynamics
            % From "Feedback Control of Dynamic Bipedal Robot Locomotion"
            theta_1 = x(1);
            theta_2 = x(2);
            theta_3 = x(3);
            dtheta_1 = x(4);
            dtheta_2 = x(5);
            dtheta_3 = x(6);
            theta = [theta_1;theta_2;theta_3];
            dtheta = [dtheta_1;dtheta_2;dtheta_3];
            q = [1 0 -1;
                 0 1 -1;
                 0 0  1] * theta + [pi;pi;0];
            qd = [1 0 -1;
                  0 1 -1;
                  0 0  1] * dtheta;
            qd1 = qd(1);
            qd2 = qd(2);
            qd3 = qd(3);
            
            
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
            C(1,2) = -0.5 * obj.m * obj.r * obj.r * sin(theta_1 - theta_2) * qd2;
            C(1,3) = obj.M_T * obj.r * obj.l * sin(theta_1 - theta_3) * qd3;
            C(2,1) = 0.5 * obj.m * obj.r * obj.r * sin(theta_1 - theta_2) * qd1;
            C(2,2) = 0.0;
            C(2,3) = 0.0;
            C(3,1) = -obj.M_T * obj.r * obj.l * sin(theta_1 - theta_3) * qd1;
            C(3,2) = 0.0;
            C(3,3) = 0.0;
        end
    end
end

