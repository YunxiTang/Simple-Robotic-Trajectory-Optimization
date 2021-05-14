%%% relaxed barrier function
clc;
clear;
x1 = -1:0.00001:20;
x2 = 0: 0.00001:20;
x3 = 0:0.1:20;
N = length(x1);
y1 = zeros(N,1);
for i=1:N
    z = x1(i);
    y1(i) = 0.2*relaxed_barrier(z, 0.1);
end

y2 = 0.2*barrier(x2);
y3 = 1 ./ x3;

%%
figure(1);
plot(x1,y1,'Color',[0.8500 0.3250 0.0980],'LineWidth',4.0);hold on;
plot(x2,y2,'Color',[0.3010 0.7450 0.9330],'LineWidth',3.0,'LineStyle','-.'); hold on;
plot(x3,y3,'Color',[0.4940 0.1840 0.5560],'LineWidth',3.0); hold on;
plot(x2,0*x2,'k--','LineWidth',2.0);hold on;
plot(0*x2,linspace(0,15,numel(x2)),'k--','LineWidth',2.0);hold on;

xlabel('$h$','Interpreter','latex','FontSize',20);
ylabel('$B(h)$','Interpreter','latex','FontSize',20);
legend('Relaxed Log Barrier Function($\mu=0.2$,$\delta=0.1$)','Log Barrier Function','$1/x$ Barrier Function','Indicator Function',...
       'Interpreter','latex','Fontsize',15);
% legend('boxoff')
axis equal;
grid off;

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
