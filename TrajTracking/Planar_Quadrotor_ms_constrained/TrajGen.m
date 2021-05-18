function [tref,xref,uref] = TrajGen()
%TRAJGEN Trajectory Generation
    tref = linspace(0.0, 6.0, 601);
    xref = zeros(6, 601);
    uref = zeros(2, 600);
    xref(1,:) = 3 .* tref;
    xref(2,:) = 2.5 + sin(2 * pi * tref);
end

