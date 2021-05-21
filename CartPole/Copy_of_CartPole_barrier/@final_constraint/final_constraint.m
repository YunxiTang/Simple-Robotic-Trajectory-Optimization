classdef final_constraint < handle
    %FINAL_CONSTRAINT final_state_constraint here
    
    properties
        delta = 0.1,
        mu = 1.1,
        t = 10,
        N_ineq = 0
    end
    
    methods
        function obj = final_constraint(Nineq, t, mu, delta)
            %FINAL_CONSTRAINT Construct an instance of this class
            fprintf('[INFO]: Creating Final Constraint for Optimization Problem. \n');
            obj.N_ineq = Nineq;
            if nargin > 1
                obj.t = t;
                obj.mu = mu;
                obj.delta = delta;
            end
        end
        
        function [] = update_t(obj)
            obj.t = min(obj.t * obj.mu, 1e5); 
            obj.delta = max(obj.delta / 1.01, 1e-4);
        end
        
        function Penalty = penalty(obj, x)
            % qudratic extension of RLB
            cons_value = obj.final_c(x);
            Nc = obj.N_ineq;
            pent = zeros(Nc, 1);
            for i=1:Nc
                h = -cons_value(i);
                if h > obj.delta
                    pent(i) = 1/obj.t * (-log(h));
                else
                    pent(i) = 1/obj.t * (1/2*(((h-2*obj.delta)/obj.delta).^2-1)-log(obj.delta));
                end
            end
            Penalty = sum(pent);
        end
        
        function Penalty = e_penalty(obj, x)
            % expential extension of RLB
            cons_value = obj.final_c(x);
            Nc = obj.N_ineq;
            pent = zeros(Nc, 1);
            for i=1:Nc
                h = -cons_value(i);
                if h > obj.delta
                    pent(i) = 1/obj.t * (-log(h));
                else
                    pent(i) = 1/obj.t * (exp(1-h/obj.delta)-1-log(obj.delta));
                end
            end
            Penalty = sum(pent);
        end
        
        function [Penal, Penal_x, Penal_xx] = penalty_info(obj, x)
            C = obj.final_c(x);
            Penal = obj.penalty(x);
            Nc = obj.N_ineq;
            Nx = numel(x);
            Penal_x = zeros(Nx, 1);
            Penal_xx = zeros(Nx, Nx);
            Cx = obj.final_cx(x);
            for i=1:Nc
                Ci = C(i);
                if -Ci > obj.delta
                    Penal_x = Penal_x + 1/obj.t * (- 1 / Ci * Cx(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * 1 / (Ci^2) * Cx(i,:)'*Cx(i,:);
                    
                else
                    T = (2*obj.delta + Ci) / (obj.delta^2);
                    Penal_x = Penal_x + 1/obj.t * (T * Cx(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * 1 / obj.delta^2 * Cx(i,:)'*Cx(i,:);
                end
            end
        end
        
        function [Penal, Penal_x, Penal_xx] = e_penalty_info(obj, x)
            C = obj.final_c(x);
            Penal = obj.e_penalty(x);
            Nc = obj.N_ineq;
            Nx = numel(x);
            Penal_x = zeros(Nx, 1);
            Penal_xx = zeros(Nx, Nx);
            Cx = obj.final_cx(x);
            for i=1:Nc
                Ci = C(i);
                if -Ci > obj.delta
                    Penal_x = Penal_x + 1/obj.t * (- 1 / Ci * Cx(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * 1 / (Ci^2) * Cx(i,:)'*Cx(i,:);    
                else
                    Penal_x = Penal_x + 1/obj.t * (exp(1-(-Ci/obj.delta)) * 1/obj.delta * Cx(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * exp(1-(-Ci/obj.delta)) * 1/obj.delta^2 * Cx(i,:)' * Cx(i,:);
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

