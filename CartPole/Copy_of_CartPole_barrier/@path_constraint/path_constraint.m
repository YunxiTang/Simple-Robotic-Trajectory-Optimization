classdef path_constraint < handle
    %path_CONSTRAINT path_constraint
    
    properties
        delta = 0.1,
        mu = 1.2,
        t = 20,
        N_ineq = 0
    end
    
    methods
        function obj = path_constraint(Nineq, t, mu, delta)
            %CONSTRAINT Construct an instance of this class
            fprintf('[INFO]: Creating Path Constraint for Optimization Problem. \n');
            obj.N_ineq = Nineq;
            if nargin > 1
                obj.t = t;
                obj.mu = mu;
                obj.delta = delta;
            end
        end
        
        function [] = update_t(obj)
            obj.t = min(obj.t * obj.mu, 1e5);  
            obj.delta = max(obj.delta / 1.1, 1e-4);
        end
        
        function Penalty = penalty(obj, x, u)
            cons_value = obj.c(x, u);
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
        
        function Penalty = e_penalty(obj, x, u)
            cons_value = obj.c(x, u);
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
        
        function [Penal, Penal_x, Penal_u, Penal_xu, Penal_ux, Penal_xx, Penal_uu] = penalty_info(obj, x, u)
            C = obj.c(x, u);
            Penal = obj.penalty(x, u);
            Nc = length(C);
            Nx = numel(x);
            Nu = numel(u);
            Penal_x = zeros(Nx, 1);
            Penal_u = zeros(Nu, 1);
            Penal_xu = zeros(Nx, Nu);
            Penal_ux = zeros(Nu, Nx);
            Penal_xx = zeros(Nx, Nx);
            Penal_uu = zeros(Nu, Nu);
            Cx = obj.cx(x, u);
            Cu = obj.cu(x, u);
            for i=1:Nc
                Ci = C(i);
                if -C(i) > obj.delta
                    Penal_x = Penal_x + 1/obj.t * (- 1 / Ci * Cx(i,:)');
                    Penal_u = Penal_u + 1/obj.t * (- 1 / Ci * Cu(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * 1 / (Ci^2) * Cx(i,:)'*Cx(i,:);
                    Penal_uu = Penal_uu + 1/obj.t * 1 / (Ci^2) * Cu(i,:)'*Cu(i,:);
                    
                else
                    T = (2*obj.delta + Ci) / (obj.delta^2);
                    Penal_x = Penal_x + 1/obj.t * (T * Cx(i,:)');
                    Penal_u = Penal_u + 1/obj.t * (T * Cu(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * 1 / obj.delta^2 * Cx(i,:)'*Cx(i,:);
                    Penal_uu = Penal_uu + 1/obj.t * 1 / obj.delta^2 * Cu(i,:)'*Cu(i,:);
                end
            end
        end
        
        function [Penal, Penal_x, Penal_u, Penal_xu, Penal_ux, Penal_xx, Penal_uu] = e_penalty_info(obj, x, u)
            C = obj.c(x, u);
            Penal = obj.e_penalty(x, u);
            Nc = length(C);
            Nx = numel(x);
            Nu = numel(u);
            Penal_x = zeros(Nx, 1);
            Penal_u = zeros(Nu, 1);
            Penal_xu = zeros(Nx, Nu);
            Penal_ux = zeros(Nu, Nx);
            Penal_xx = zeros(Nx, Nx);
            Penal_uu = zeros(Nu, Nu);
            Cx = obj.cx(x, u);
            Cu = obj.cu(x, u);
            for i=1:Nc
                Ci = C(i);
                if -C(i) > obj.delta
                    Penal_x = Penal_x + 1/obj.t * (- 1 / Ci * Cx(i,:)');
                    Penal_u = Penal_u + 1/obj.t * (- 1 / Ci * Cu(i,:)');
                    Penal_xx = Penal_xx + 1/obj.t * 1 / (Ci^2) * Cx(i,:)'*Cx(i,:);
                    Penal_uu = Penal_uu + 1/obj.t * 1 / (Ci^2) * Cu(i,:)'*Cu(i,:);
                    
                else
                    Penal_x = Penal_x + 1/obj.t * (exp(1-(-Ci/obj.delta)) * 1/obj.delta) * Cx(i,:)';
                    Penal_xx = Penal_xx + 1/obj.t * (exp(1-(-Ci/obj.delta)) * 1/obj.delta^2) * Cx(i,:)' * Cx(i,:);
                    Penal_u = Penal_u + 1/obj.t * (exp(1-(-Ci/obj.delta)) * 1/obj.delta) * Cu(i,:)';
                    Penal_uu = Penal_uu + 1/obj.t * (exp(1-(-Ci/obj.delta)) * 1/obj.delta^2) * Cu(i,:)' * Cu(i,:);
                end
            end
        end
        
    end
    
    methods (Static)
        % compute constraint
        c = c(in1,in2);
        
        % compute gradients
        cu = cu(in1,in2);
        cx = cx(in1,in2)
    end
end

