classdef newsosbalance < sosbalance
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stop_critia
    end
    
    methods
        function newsosbalance = newsosbalance(in1, varargin)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            % assign the superclass portion
            newsosbalance = newsosbalance@sosbalance(varargin{:});
            if nargin > 0
                newsosbalance.stop_critia = in1;
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

