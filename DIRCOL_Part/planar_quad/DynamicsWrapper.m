function [dx] = DynamicsWrapper(t,x,u,dyn_func)
%DYNAMICSWRAPPER dynamics warpper for TrajOpt

    [Nx, Nt] = size(x);
    dx = zeros(Nx, Nt);
    for i=1:Nt
        dx(:,i) = dyn_func(t(i), x(:,i), u(:,i));
    end
end

