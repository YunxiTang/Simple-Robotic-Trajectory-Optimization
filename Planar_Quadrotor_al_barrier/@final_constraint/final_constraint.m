classdef final_constraint < handle
    %FINAL_CONSTRAINT final_state_constraint here
    
    properties
        delta = 0.1,
        mu = 100,
        t = 1
    end
    
    methods
        function obj = final_constraint(t, mu, delta)
            %FINAL_CONSTRAINT Construct an instance of this class
            fprintf('[INFO]: Creating Final Constraint for Optimization Problem. \n');
            if nargin > 1
                obj.t = t;
                obj.mu = mu;
                obj.delta = delta;
            end
        end
        
        function [] = update_t(obj)
            obj.t = min(obj.t + obj.mu, 1e5);
            obj.delta  = max(obj.delta / 1.001, 1e-9);         
        end
        
        function Penalty = penalty(obj, x)
            cons_value = obj.final_c(x);
            Nc = length(cons_value);
            pent = zeros(Nc, 1);
            for i=1:Nc
                if -cons_value(i) > obj.delta
                    pent(i) = 1/obj.t * (-log(-cons_value(i)));
                else
                    h = cons_value(i);
                    pent(i) = 1/obj.t * (1/2*(((-h-2*obj.delta)/obj.delta).^2-1)-log(obj.delta));
                end
            end
            Penalty = sum(pent);
        end
        
        function [Penal, Penal_x, Penal_xx] = penalty_info(obj, x)
            C = obj.final_c(x);
            Penal = obj.penalty(x);
            Nc = length(C);
            Nx = numel(x);
            Penal_x = zeros(Nx, 1);
            Penal_xx = zeros(Nx, Nx);
            Cx = obj.final_cx(x);
            for i=1:Nc
                Ci = C(i);
                if -C(i) > obj.delta
                    Penal_x = Penal_x + 1/obj.t * (- 1 / Ci * Cx(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * 1 / (Ci^2) * Cx(i,:)'*Cx(i,:);
                    
                else
                    T = (2*obj.delta + Ci) / (obj.delta^2);
                    Penal_x = Penal_x + 1/obj.t * (T * Cx(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * 1 / obj.delta^2 * Cx(i,:)'*Cx(i,:);
                end
            end
        end
    end
    
    methods (Static)
        % compute constraint
        c = final_c(in1);
        
        % compute gradients
        cx = final_cx(in1)
    end
end

