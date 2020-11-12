%%% consider minmize an non-convex polynominal
x = sdpvar(1,1);
y = sdpvar(1,1);
rho = sdpvar(1,1);

p1 = 4*x^2 - 21/10*x^4 + 1/3*x^6 + x*y -4*y^2 + 4*y^4;
F = [sos(p1-rho)];
[sol,m,B,residual] = solvesos(F,-rho,[],[rho]);

%% visualization
x = linspace(-2.5,2.5,500);
y = linspace(-2.5,2.5,500);
[X,Y]=meshgrid(x,y);
P = 4*X.^2 - 21/10*X.^4 + 1/3*X.^6 + X.*Y -4*Y.^2 + 4*Y.^4;
[m,c]= contour(X,Y,P,100);
c.LineWidth = 2.0;
axis equal;