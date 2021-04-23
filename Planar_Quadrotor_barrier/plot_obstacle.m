function Flag = plot_obstacle(obstacle, map)
if nargin > 1
    figure(map);
end
aplha=0:0.01:2*pi;
n_ineq = size(obstacle, 2);
for i=1:n_ineq
    % the i-th circle
    r = obstacle(3,i); 
    xi_c = obstacle(1,i);
    yi_c = obstacle(2,i);
    cx = xi_c + r*cos(aplha);
    cy = yi_c + r*sin(aplha);
    clr = abs([sin(i/1.2) sin(i/1.5) sin(i)]);
    plot(cx, cy, '-','Color',clr); hold on;
    h = fill(cx, cy, clr);
    set(h,'edgealpha',0.1,'facealpha',0.7);
    axis equal;
end
xlabel('$x$','Interpreter','latex','FontSize',15);
ylabel('$y$','Interpreter','latex','FontSize',15);
title("$2D\;Map$", 'Interpreter','latex','FontSize',15);
grid on;
Flag = 1;
end
