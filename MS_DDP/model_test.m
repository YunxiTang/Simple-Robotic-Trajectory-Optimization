%%% TEST DYNAMICS
clc;clear;
pendubot = Pendubot();

%% simulation test
dt = 0.001;
t_end = 30;
t_span = [0.0, t_end];
x0 = [pi/4;0.0;0.0;0.0];
x = x0;
X = [];
T = [];
for t = t_span(1):dt:t_span(2)
    x_s = pendubot.rk45(x,[0], dt);
    x = x_s;
    X = [X x_s];
    T = [T t];
end
u = @(t,x)(0.0);

opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
sol = ode45(@(t,x)pendubot.Dynamics(t,x,u(t,x)),[0:dt:t_end],x0,opts);
t = sol.x;
state = sol.y;
%%
figure(111);
for i=1:4
    subplot(4,1,i);
    plot(T,X(i,:),'g','LineWidth',3.5);hold on; plot(t, state(i,:),'k-.','LineWidth',1.5);
    grid on;
    legend('$RK45$','$ODE45$','Interpreter','latex','FontSize',11);
    tit = ['$x_', num2str(i),'$'];
    title(tit,"Interpreter", "latex", "FontSize", 15);
end
xlabel('$Time [s]$','Interpreter','latex','FontSize',12);
hold off;
