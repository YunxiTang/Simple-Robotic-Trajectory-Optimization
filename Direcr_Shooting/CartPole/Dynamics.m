function [dx] = Dynamics(x,u)
%DYNAMICS 
    % dynamics
    % model params
    m1 = 2.0;
    m2 = 0.5;
    l = 0.5;
    g = 9.81;

    q = x(1);
    theta = x(2);
    dq = x(3);
    dtheta = x(4);
    
    M = [cos(theta)         l;
         m1 + m2    m2*l*cos(theta)];
    G = [-g * sin(theta);
         u + m2 * l * dtheta^2 * sin(theta)];
    acc = M \ G;
    dx = [dq;
          dtheta;
          acc]; 
end

