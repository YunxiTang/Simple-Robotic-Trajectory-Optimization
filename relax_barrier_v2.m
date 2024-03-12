%%% relaxed barrier function
clc;
clear;
ub = 10;
x1 = -0.8:0.00001:(ub+0.8);
x2 = 0: 0.00001:ub;
x3 = 0:0.1:ub;
N = length(x1);
y1 = zeros(N,1);
for i=1:N
    z = x1(i);
    y1(i) = 0.2*relaxed_barrier(z, 0.1) + 0.2*relaxed_barrier(ub-z, 0.1);
end

y2 = 0.2*barrier(x2) + 0.2*barrier(ub-x2);
y3 = 1 ./ x3 + 1 ./ (ub - x3);

%%
figure(1);

plot(x1,y1,'r','LineWidth',3.0);hold on;

plot(x2,y2,'b','LineWidth',3.0,'LineStyle','-.'); hold on;

plot(x3,y3,'g','LineWidth',3.0); hold on;

plot(x2,0*x2,'k-.','LineWidth',2.0);hold on;
plot(0*x2,linspace(0,10,numel(x2)),'k-.','LineWidth',2.0);hold on;
plot(0*x2+ub,linspace(0,10,numel(x2)),'k-.','LineWidth',2.0);hold on;
xlabel('$c$','Interpreter','latex','FontSize',20);
ylabel('$B(c)$','Interpreter','latex','FontSize',20);
legend('Relaxed Log Barrier Function($\mu=0.2$,$\delta=0.1$)','Log Barrier Function','$1/x$ Barrier Function','Indicator Function',...
       'Interpreter','latex','Fontsize',15);
% legend('boxoff')
axis equal;
grid on;

%%
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
