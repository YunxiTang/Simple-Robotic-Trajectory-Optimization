c1 = [2.2547
   -0.7851
    1.4554];
x1 = -2.5:0.01:2.5;
x2 = x1;
[X1, X2] = meshgrid(x1, x2);
V = X1.^2*c1(1)+X1.*X2*c1(2)+X2.^2*c1(3);
surf(X1,X2,V)
condition1 = V <= 0.7;
condition2 = V >= 0;
output = ones(length(x1),length(x2));
output(~(condition1&condition2)) = 0;
imshow(output,'xdata',x1,'ydata',x2);
axis equal;
axis on;
xlabel('x1');
ylabel('x2');