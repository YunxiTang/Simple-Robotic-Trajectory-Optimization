classdef CartPole
    %CARTPOLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m_c,    % mass of cart
        m_p,    % mass of pole
        m_l,    % length of pole
    end
    
    methods
        function obj = CartPole(inputArg1,inputArg2)
            %CARTPOLE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

