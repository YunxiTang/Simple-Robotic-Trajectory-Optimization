%%%  Van Der Pol Oscillator: Lyapunov function search iteratively
clear all;
clc;
mu = 1.0;
r = 2.0;
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x = [x1;x2];

f = [-x2;
     -mu*(1-x1^2)*x2+x1];
 
[V, c1] = polynomial([x1;x2],2,2);

[L, c2] = polynomial([x1;x2],2);

norm22 = 1e-7*(x1^2+x2^2);

Expr = -jacobian(V,[x1;x2])*f - L*(r - V);

F = [sos(V), sos(L), sos(Expr)];

h=solvesos(F,[],[],[c1;c2]);

