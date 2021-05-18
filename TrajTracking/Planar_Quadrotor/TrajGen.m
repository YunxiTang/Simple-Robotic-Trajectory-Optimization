function [tref,xref,uref] = TrajGen()
%TRAJGEN Trajectory Generation
    tref = linspace(0.0, 10.0, 1001);
    xref = zeros(6, 1001);
    uref = zeros(2, 1000);
    xref(1,:) = 0.4 .* tref;
    xref(2,:) = min(max(4 + 3 * sin(0.2 * pi * tref), 3.5),4.5);
    xref(4,:) = 0.4 * ones(1001, 1);
    for i=1:1001
        if xref(2,i) == 3.5 || xref(2,i) == 4.5
            xref(5,i) = 0.0;
        else
            xref(5,i) = 3 * 0.2 * pi * cos(0.2 * pi * tref(i));
        end
    end
end

