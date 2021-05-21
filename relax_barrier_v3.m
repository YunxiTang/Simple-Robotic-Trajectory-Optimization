%%% relaxed barrier function
clc;
clear;
x1 = -1:0.00001:2;
x2 = 0: 0.00001:2;
% x3 = 0:0.1:20;
N = length(x1);
y1 = zeros(N,1);
y2 = zeros(N,1);
y3 = zeros(N,1);
y4 = zeros(N,1);
y5 = zeros(N,1);
for i=1:N
    z = x1(i);
    y1(i) = 0.1*relaxed_barrier(z, 0.1);
    y2(i) = 0.01*relaxed_barrier(z, 0.01);
    y3(i) = 0.001*relaxed_barrier(z, 0.001);
    y4(i) = 1e-4*relaxed_barrier(z, 1e-4);
    y5(i) = 1e-5*relaxed_barrier(z, 1e-5);
end

%%
figure(1);
plot(x1,y1,'c-.','LineWidth',2.4);hold on;
plot(x1,y2,'m-.','LineWidth',2.4);hold on;
plot(x1,y3,'g-.','LineWidth',2.4);hold on;
plot(x1,y4,'b-.','LineWidth',2.4);hold on;
plot(x1,y5,'r-.','LineWidth',2.4);hold on;
plot(x2,0*x2,'k--','LineWidth',2.0);hold on;
plot(0*(0:0.1:2),0:0.1:2,'k--','LineWidth',2.0);hold on;

xlabel('$c$','Interpreter','latex','FontSize',20);
ylabel('$B(c)$','Interpreter','latex','FontSize',20);
legend('$\mu=0.1$, $\delta=0.1$',...
       '$\mu=0.01$, $\delta=0.01$',...
       '$\mu=0.001$, $\delta=0.001$', ...
       '$\mu=0.0001$, $\delta=0.0001$',...
       '$\mu=0.00001$, $\delta=0.00001$',...
       'Indicator Function',...
       'Interpreter','latex','Fontsize',15);
title('Relaxed Log Barrier Functions','Interpreter','latex','Fontsize',20);
legend('boxoff')
xlim([-0.1 1]);
ylim([0 0.5]);
% axis equal;
grid off;

%%
function value = relaxed_barrier(x, delta)
    if delta < x
        value = -log(x);
    else
%         value = 1/2*(((x-2*delta)/delta).^2-1)-log(delta);
        value = exp(1-x/delta)-1-log(delta);
    end  
end