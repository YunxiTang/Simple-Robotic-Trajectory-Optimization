function [f] = VanderDyna(t,x)
%VANDERDYNA Summary of this function goes here
    x1 = x(1);
    x2 = x(2);
    f = [-x2;
         x1-x2+x2^3];
end

