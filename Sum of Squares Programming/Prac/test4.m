%% Motzkin Polynomial
x = -1.15:0.01:1.15;
y = x;

[X, Y] = meshgrid(x, y);
M = X.^2.*Y.^4+X.^4.*Y.^2 + 1 - 3*X.^2.*Y.^2;
surf(X,Y,M);