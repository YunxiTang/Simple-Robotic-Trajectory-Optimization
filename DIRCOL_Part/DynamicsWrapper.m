function [dx] = DynamicsWrapper(dyn_func,x,u)
%DYNAMICSWRAPPER dynamics warpper for TrajOpt
    [Nx, Nt] = size(x);
    dq = x((Nx/2+1):Nx,:);
    ddq = zeros(Nx/2, Nt);
    for i=1:Nt
        ddq(:,i) = dyn_func(0.0, x(:,i), u(:,i));
    end
    dx = [dq;
          ddq];
end

