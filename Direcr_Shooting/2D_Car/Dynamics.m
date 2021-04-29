function [dq] = Dynamics(q,u)
%DYNAMICS 
    % dynamics

    x = q(1);
    y = q(2);
    theta = q(3);
    v = q(4);
    
    u_theta = u(1);
    u_v = u(2);
    
    dq = [v*sin(theta);
          v*cos(theta);
          u_theta * v;
          u_v]; 
end

