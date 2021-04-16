classdef Bipedal < handle
    %BIPEDAL Three Link Bipedal Robot
    
    properties
        Name = 'Bipedal',
        l = 0.3,      % torso length
        r = 0.4,      % leg length
        M_T = 0.3,    % torso mass
        M_H = 0.4,    % hip mass
        m = 0.2,      % leg mass
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
            M(1,1) = (5 * obj.m / 4 + obj.M_H + obj.M_T);
        end
    end
end

