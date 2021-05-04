clc;clear;
M1 = load('M1.mat');        S1 = M1.J_hist;     L1 = length(S1);
M2 = load('M2.mat');        S2 = M2.J_hist;     L2 = length(S2);
M5 = load('M5.mat');        S5 = M5.J_hist;     L5 = length(S5);
M10 = load('M10.mat');      S10 = M10.J_hist;   L10 = length(S10);
M20 = load('M20.mat');      S20 = M20.J_hist;   L20 = length(S20);
M25 = load('M25.mat');      S25 = M25.J_hist;   L25 = length(S25);
M50 = load('M50.mat');      S50 = M50.J_hist;   L50 = length(S50);
M100 = load('M100.mat');    S100 = M100.J_hist; L100 = length(S100);
M250 = load('M250.mat');    S250 = M250.J_hist; L250 = length(S250);
M500 = load('M500.mat');    S500 = M500.J_hist; L500 = length(S500);
t1 = load('T1.mat');        te1 = t1.telapsed;
t2 = load('T2.mat');        te2 = t2.telapsed;
t5 = load('T5.mat');        te5 = t5.telapsed;
t10 = load('T10.mat');      te10 = t10.telapsed;
t20 = load('T20.mat');      te20 = t20.telapsed;
t25 = load('T25.mat');      te25 = t25.telapsed;
t50 = load('T50.mat');      te50 = t50.telapsed;
t100 = load('T100.mat');    te100 = t100.telapsed;
t250 = load('T250.mat');    te250 = t250.telapsed;
t500 = load('T500.mat');    te500 = t500.telapsed;

figure(111);
points = 1;
plot(1:points:L1,S1(1:points:L1),'r-o','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L2,S2(1:points:L2),'b-o','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L5,S5(1:points:L5),'g-o','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L10,S10(1:points:L10),'b-x','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L20,S20(1:points:L20),'r-x','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L25,S25(1:points:L25),'g-x','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L50,S50(1:points:L50),'m-x','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L100,S100(1:points:L100),'m-d','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L250,S250(1:points:L250),'r-d','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(1:points:L500,S500(1:points:L500),'b-d','LineWidth',1.0,'MarkerSize',8.0); hold on;
plot(S1(end)*ones(length(S1)),'k--','LineWidth',2.0);
ha=gca;
set(ha,'yscale','log');
set(ha,'xscale','log');
grid on;
title('Cost v.s. Iterations','Interpreter','latex','FontSize',15);
xlabel('Iter','Interpreter','latex','FontSize',15);
ylabel('Cost','Interpreter','latex','FontSize',15);
lgd = legend('Shooting Phase: 1','Shooting Phase: 2','Shooting Phase: 5','Shooting Phase: 10','Shooting Phase: 20','Shooting Phase: 25','Shooting Phase: 50','Shooting Phase: 100',...
       'Shooting Phase: 250','Shooting Phase: 500','Interpreter','latex','FontSize',12);
lgd.NumColumns = 2;
   
fig = figure(222);
Ms = [1 2 5 10 20 25 50 100 250 500];
Iters = [L1 L2 L5 L10 L20 L25 L50 L100 L250 L500];
time = [te1 te2 te5 te10 te20 te25 te50 te100 te250 te500];

set(fig,'defaultAxesColorOrder',[[1 0 0]; [0 0 0]]);
yyaxis left
semilogx(Ms,Iters,'rd--','LineWidth',0.5,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8.0);hold on;
ylabel('Iter','Interpreter','latex','FontSize',15);
yyaxis right
semilogx(Ms,time, 'ks--','LineWidth',0.5,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8.0);
ylabel('Computation Time [s]','Interpreter','latex','FontSize',15);
title('Different  Shooting Phase','Interpreter','latex','FontSize',15);
xlabel('Shooting Phase','Interpreter','latex','FontSize',15);

grid on;