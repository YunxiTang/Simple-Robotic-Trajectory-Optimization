%%% Van Der Pol Oscillator (iteratively)
clear all;clc;
syms x y;
vars = [x;y];
prog = sosprogram(vars);
Z = monomials(vars,0:4);
% f = dx/dt
f = [-x-(x^2-1)*y;
     y];

% V
[prog, V] = sospolyvar(prog, Z, 'wscoeff');

% constraint 1: V is SOS
prog = sosineq(prog, V-x^2-y^2);

dV = diff(V,x)*f(1)+diff(V,y)*f(2)-x^2-y^2;

% constraint 2:-dV is SOS in B(1)
prog = sosineq(prog,-dV);

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
Sol_V = sosgetsol(prog,V);
