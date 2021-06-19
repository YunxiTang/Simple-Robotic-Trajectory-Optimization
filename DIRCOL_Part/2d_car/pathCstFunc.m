function [Cp_ineq,Cp_eq] = pathCstFunc(t, x, u, params)
%PATHCSTFUNC path constraint function
%   Cp_ineq: path constraint / inequality constraint (<=0)
%   Cp_eq:   path constraint /   equality constraint (==0)
    nTimes = size(t,2);
    obstacle = params.Obstacles;
    nIneq = size(obstacle,1);
    Cp_ineq = zeros(nIneq,nTimes);
    % no equality constraints
    Cp_eq = [];
    
    for i=1:nTimes
        xi = x(:,i);
        for j=1:nIneq
            Cp_ineq(j,i) = obstacle(3,j)*obstacle(3,j) - ...
                           ((xi(1)-obstacle(1,j))*(xi(1)-obstacle(1,j)) + ...
                            (xi(2)-obstacle(2,j))*(xi(2)-obstacle(2,j)));
        end
    end
end

