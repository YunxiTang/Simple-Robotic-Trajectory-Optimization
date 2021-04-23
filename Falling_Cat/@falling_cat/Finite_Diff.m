function [Jx, Ju] = Finite_Diff(func, qbar, ubar, dt, interval)
    if nargin == 5
        h = interval;
    else
        h = 1e-3;
    end
    Nx = numel(qbar);
    Nu = numel(ubar);
    Jx = zeros(Nx,Nx);
    Ju = zeros(Nu,Nu);
    Hx = eye(Nx) * h;
    Hu = eye(Nu) * h;
    for i=1:Nx
        Jx(:,i) = (func(qbar + Hx(:,i), ubar, dt) - func(qbar - Hx(:,i), ubar, dt)) / (2*h);
    end
    
    for i=1:Nu
        Ju(:,i) = (func(qbar, ubar+ Hu(:,i), dt) - func(qbar , ubar - Hu(:,i), dt)) / (2*h);
    end
end
