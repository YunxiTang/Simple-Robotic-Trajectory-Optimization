classdef sosbalance
    % sos balance controller
    % public and protected: access control
    properties(Access = public)
        Lyafunc
        controller
        LagMul
    end
    
    properties(Access = protected)
        x = sdpvar(1,1);
        y = sdpvar(1,1);
   
    end
    
    % constant and dependent: modification control
    properties(Constant)
        epsilon = 1e-6;
    end
    
    properties(Dependent)
        state
    end
    
    methods
        function obj = sosbalance(Lya,controller,L)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 3
                obj.Lyafunc = Lya;
                obj.controller = controller;
                obj.LagMul = L;
            end
        end
        
        function [] = show(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            disp('Current Lyapunov funtion is');
            disp(obj.Lyafunc);
        end
        
        function state = get.state(obj)
            % dependent function (whenever we call obj.state, run this function directly)
            state = obj.x^2 + obj.y^2;
        end
    end
end

