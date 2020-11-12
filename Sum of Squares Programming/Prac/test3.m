%%% YALMIP and SOSTOOL to test p(x)is SOS or not
clear all;clc;
syms x y;
vars = [x;y];
p = 4*x^4 + 4*x^3*y -7*x^2*y^2 - 2*x*y^3 + 10*y^4;
prog = sosprogram(vars);
prog = sosineq(prog,p);
[h,INFO] = sossolve(prog);

%% YALMIP
clear all;clc;
sdpvar x y;
p = 4*x^4 + 4*x^3*y -7*x^2*y^2 - 2*x*y^3 + 10*y^4;
F = [];
F = [F;sos(p)];
solvesos(F);
sdisplay(p);