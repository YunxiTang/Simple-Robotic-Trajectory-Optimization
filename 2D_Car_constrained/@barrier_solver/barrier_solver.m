classdef barrier_solver < handle
    %barrier solver to polish the coarse solution from AL-DDP
    
    properties
        Property1
    end
    
    methods
        function obj = barrier_solver(inputArg1,inputArg2)
            %QP_SOLVER Construct an instance of this class
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

