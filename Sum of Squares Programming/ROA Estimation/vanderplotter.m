q1 = linspace(-1,1,100);
q2 = q1;
[Q1,Q2] = meshgrid(q1,q2);
Lya = 1.5*Q1.^2+Q2.^2-Q1.*Q2;
h = value(gamma);

figure;
surf(Q1,Q2,Lya);hold on;
fill3([-1;-1;1;1],[-1;1;1;-1],[h;h;h;h],'r');
xlabel('q1');
ylabel('q2');
zlabel('Lypunov V_{ini}');
axis equal;
figure;
contour(Q1,Q2,Lya);
