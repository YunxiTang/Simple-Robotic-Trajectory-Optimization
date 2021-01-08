%%% TEST DYNAMICS
clc;clear;
param = struct('l1',1.0,'l2',1.0,...
               'm1',1.0,'m2',1.0,...
               'lc1',1.5,'lc2',0.5,...
               'b1',0.1,'b2',0.1,...
               'I1',0.33,'I2',0.33,...
               'g',9.81);

acrobot = Acrobot_mdl(param);

%% simulation test
dt = 0.001;
t_span = [0.0, 10.0];
x0 = [pi;0.0;0.0;0.0];
x = x0;
X = [];
T = [];
for t = t_span(1):dt:t_span(2)
    x_s = acrobot.rk45(x, [0], dt);
    x = x_s;
    X = [X x_s];
    T = [T t];
end
figure(111);
for i=1:4
    subplot(4,1,i);
    plot(T,X(i,:),'LineWidth',2.0);
    tit = ['$x_', num2str(i),'$'];
    title(tit,"Interpreter", "latex", "FontSize", 20);
end
hold off;
acrobot.run_animation(T,X);
