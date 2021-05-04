function [xd] = Dynamics(x,u)
%DYNAMICS 
    % dynamics
    m = 0.8;
    g = 9.81;
    l = 0.3;
    J = 1.0;
    
    theta = x(3);
    ddx = -(1 / m) * (u(1) + u(2)) * sin(theta);
    ddy =  (1 / m) * (u(1) + u(2)) * cos(theta) - g;
    ddtheta =  (1 / J) * (l / 2) * (u(2) - u(1));
    xd = [x(4);
          x(5);
          x(6);
          ddx;
          ddy;
          ddtheta];
end

