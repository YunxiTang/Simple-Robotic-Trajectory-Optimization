function [Cp_ineq,Cp_eq] = pathCstFunc(t, x, u, params)
%PATHCSTFUNC path constraint function
%   Cp_ineq: path constraint / inequality constraint (<=0)
%   Cp_eq:   path constraint /   equality constraint (==0)
    nTimes = size(t,2);
    nIneq = params.NJ;
    Cp_ineq = zeros(nIneq,nTimes);
    % no equality constraints
    Cp_eq = [];
    
    for i=1:nTimes
        xi = x(:,i);
        for j=1:nIneq
            Cp_ineq(j,i) = xi(1) - pi;
        end
    end
end

