%%% Lypunov search for given stable system (global stable)
clear all;
clc;
syms x y;
vars = [x;y];
prog = sosprogram(vars);

% dx/dt
f = [-y-3/2*x^2-1/2*x^3;
     3*x-y];
 
Z = monomials(vars,[1,2,3,4]);

[prog,V] = sospolyvar(prog,Z,'wscoeff');

% constraint1: V is SOS
prog = sosineq(prog,V);

dV = diff(V,x)*f(1) + diff(V,y)*f(2) + x^2 + y^2;
% constraint1: -dV is SOS
prog = sosineq(prog,-dV);

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);
SOL_V= sosgetsol(prog,V);

%%
x1 = -10:0.01:10;
x2 = -10:0.01:10;
[X1,X2] = meshgrid(x1,x2);
Z = 0.1803*X1.^4 - 0.02599*X1.^3.*X2 - 0.02701*X1.^3 + 0.1776*X1.^2.*X2.^2 + 0.9783*X1.^2.*X2 + 1.695*X1.^2 + 2.182e-5*X1.*X2.^3 - 0.1272*X1.*X2.^2 - 0.5211*X1.*X2 + 1.881e-6*X1 + 0.03325*X2.^4 + 0.2298*X2.^3 + 0.6205*X2.^2 + 2.445e-6*X2;

% initial state
q01 = [5;5];q02 = [-6;6];q03 = [7;-7];q04 = [-8;-8];q05 = [-2;-2];
[t1,q1] = ode45(@(t,state)dynamics(t,state),[0 50], q01);
[t2,q2] = ode45(@(t,state)dynamics(t,state),[0 50], q02);
[t3,q3] = ode45(@(t,state)dynamics(t,state),[0 50], q03);
[t4,q4] = ode45(@(t,state)dynamics(t,state),[0 50], q04);
[t5,q5] = ode45(@(t,state)dynamics(t,state),[0 50], q05);
%%
figure(1);
[C,h]=contour(X1,X2,Z,50);
hold on;
plot(q1(:,1),q1(:,2),'linewidth',2.0);
hold on;
plot(q2(:,1),q2(:,2),'linewidth',2.0);
hold on;
plot(q3(:,1),q3(:,2),'linewidth',2.0);
hold on;
plot(q4(:,1),q4(:,2),'linewidth',2.0);
hold on;
plot(q5(:,1),q5(:,2),'linewidth',2.0);
grid on
axis equal