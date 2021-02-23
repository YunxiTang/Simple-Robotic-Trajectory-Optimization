classdef msddp_solver < base_solver
    %MSDDPSOLVER multiple shooting DDP solver
    
    properties
        M,   % segments
    end
    
    methods
        function obj = msddp_solver(inputArg1,inputArg2)
            %MSDDPSOLVER Construct an instance of this class
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

