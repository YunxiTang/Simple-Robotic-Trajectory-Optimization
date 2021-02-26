clc;clear;
M1 = load('M1.mat');        S1 = M1.J_hist;     L1 = length(S1);
M2 = load('M2.mat');        S2 = M2.J_hist;     L2 = length(S2);
M5 = load('M5.mat');        S5 = M5.J_hist;     L5 = length(S5);
M10 = load('M10.mat');      S10 = M10.J_hist;   L10 = length(S10);
M20 = load('M20.mat');      S20 = M20.J_hist;   L20 = length(S20);
M100 = load('M100.mat');    S100 = M100.J_hist; L100 = length(S100);
M250 = load('M250.mat');    S250 = M250.J_hist; L250 = length(S250);
M500 = load('M500.mat');    S500 = M500.J_hist; L500 = length(S500);
figure(111);
plot(1:5:L1,S1(1:5:L1),'r-o','LineWidth',1.0,'MarkerSize',5.0); hold on;
plot(1:5:L2,S2(1:5:L2),'b-o','LineWidth',1.0,'MarkerSize',5.0); hold on;
plot(1:5:L5,S5(1:5:L5),'g-o','LineWidth',1.0,'MarkerSize',5.0); hold on;
plot(1:5:L10,S10(1:5:L10),'b-d','LineWidth',1.0,'MarkerSize',5.0); hold on;
plot(1:5:L20,S20(1:5:L20),'r-d','LineWidth',1.0,'MarkerSize',5.0); hold on;
plot(1:5:L100,S100(1:5:L100),'m-d','LineWidth',1.0,'MarkerSize',5.0); hold on;
plot(1:5:L250,S250(1:5:L250),'r-s','LineWidth',1.0,'MarkerSize',5.0); hold on;
plot(1:5:L500,S500(1:5:L500),'b-s','LineWidth',1.0,'MarkerSize',5.0); hold on;
plot(S1(end)*ones(length(S1)),'k--','LineWidth',2.0);
ylim([15 40]);
grid on;
title('Cost v.s. Iterations','Interpreter','latex','FontSize',15);
xlabel('Iter','Interpreter','latex','FontSize',15);
ylabel('Cost','Interpreter','latex','FontSize',15);
legend('Shooting Phase: 1','Shooting Phase: 2','Shooting Phase: 5','Shooting Phase: 10','Shooting Phase: 20','Shooting Phase: 100',...
       'Shooting Phase: 250','Shooting Phase: 500','Interpreter','latex','FontSize',15);
fig = figure(222);
Ms = [1 2 5 10 20 100 250 500];
Iters = [L1 L2 L5 L10 L20 L100 L250 L500];
time = [6.802813 7.325260 6.517115 3.743957 2.845905 2.724406 3.093054 3.685360];

set(fig,'defaultAxesColorOrder',[[1 0 0]; [0 0 0]]);
yyaxis left
plot(Ms,Iters,'ro--','LineWidth',2.0,'MarkerFaceColor','r','MarkerEdgeColor','r');hold on;
ylabel('Iter','Interpreter','latex','FontSize',15);
yyaxis right
plot(Ms,time, 'ks--','LineWidth',2.0,'MarkerFaceColor','k','MarkerEdgeColor','k');
ylabel('Computation Time [s]','Interpreter','latex','FontSize',15);
title('Different  Shooting Phase','Interpreter','latex','FontSize',15);
xlabel('Shooting Phase','Interpreter','latex','FontSize',15);

grid on;