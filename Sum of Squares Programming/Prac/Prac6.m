%%% Van der Pol equations
clear all;clc;
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
epsilon = 1e-10;
rho = 0.01;

% dynamics
f = [-x2;
      x1 + (x1^2 - 1)*x2];
norm22 = epsilon * (x1^2 + x2^2);
% parameterized Lyapunov function V
[V,c] = polynomial([x1 x2],2);
[h,s] = polynomial([x1 x2],2);
dV = jacobian(V,[x1;x2])*f;
F = [sos(V-norm22),sos(-dV - h*(rho - V) - norm22)];
solvesos(F,[],[],[c;s])
