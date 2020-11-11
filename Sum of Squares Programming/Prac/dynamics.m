function [dq] = dynamics(t,q)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    x = q(1,:);
    y = q(2,:);
    dq = [-y-3/2*x^2-1/2*x^3;
           3*x-y];
end

