classdef constraint
    %CONSTRAINT build the constraints
    
    properties
        num_ineq,
        num_eq
    end
    
    methods
        function obj = constraint()
            %CONSTRAINT Construct an instance of this class
            disp('[INFO]: Creating A Template for Constraints.');
        end
        
        function outputArg = add_eqcons(obj,inputArg)
            % add equality constraints
            outputArg = obj.Property1 + inputArg;
        end
        
        function output = add_ineqcons(obj,inputArg)
        % add equality constraints
        end
        
        function c = eval_cons(obj, x, u)
        end
    end
end

