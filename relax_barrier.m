%%% relaxed barrier function
clc;
clear;
x1 = -1:0.1:1;
x2 = 0:0.00001:1;
x3 = 0:0.1:1;
N = length(x1);
y1 = zeros(N,1);
for i=1:N
    z = x1(i);
    y1(i) = 30.0*relaxed_barrier(z, 0.05);
end

y2 = barrier(x2);
y3 = 1 ./ x3;
figure(1);hold on;
plot(x1,y1,'b','LineWidth',3.0);hold on;
plot(x2,y2,'r-.','LineWidth',2.0); hold on;
plot(x3,y3,'k','LineWidth',2.0);
xlabel('x','Interpreter','latex','FontSize',14);
legend('Relaxed Log Barrier Function','Log Barrier Function','$1/x$ Barrier Function',...
       'Interpreter','latex','Fontsize',14);
grid on;


function value = relaxed_barrier(x, delta)
    if delta < x
        value = -log(x);
    else
        value = 1/2*(((x-2*delta)/delta).^2-1)-log(delta);
    end  
end

function value = barrier(x)
        value = -log(x);      
end
