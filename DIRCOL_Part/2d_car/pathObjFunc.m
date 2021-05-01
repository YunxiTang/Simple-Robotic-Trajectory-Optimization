function [J] = pathObjFunc(t,x,u,params)
%PATHOBJFUNC 
% path cost function
    nTimes = size(t,2);
    J = zeros(1,nTimes);
    for i = 1:nTimes
        ui = u(:,i);
        xi = x(:,i);
        J(i) = (ui'*params.R*ui + (xi - params.xf)' * params.Q * (xi - params.xf)) / 2;
    end
end

