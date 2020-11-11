%%% estimation of region of acttraction for 1 dimesional system
clear all; clc; close all;
syms x

vars = x;

% f = dx / dt
f = -x + x^3 ;

% Lyapunov function V
V = x^2;

prog = sosprogram(vars);

[prog,lambda] = sospolyvar(prog,[x^0;x; x^2],'wscoeff');

% constraint 1: lambda(x) SOS
prog = sosineq(prog, lambda);

% constraint 2: -dV/dt - lambda*(1-V) SOS
expr = -diff(V,x)*f - lambda*(1-V); %%it is bilinear in decision variables.
prog = sosineq(prog, expr);

solver_opt.solver = 'sedumi';

prog = sossolve(prog,solver_opt);

SOLLam = sosgetsol(prog,lambda);


