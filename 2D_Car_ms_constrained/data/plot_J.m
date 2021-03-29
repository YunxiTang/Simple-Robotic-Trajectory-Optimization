clc;clear;
M1 = load('M_1.mat');        S1 = M1.J_hist;     L1 = length(S1);
M2 = load('M_2.mat');        S2 = M2.J_hist;     L2 = length(S2);
M5 = load('M_5.mat');        S5 = M5.J_hist;     L5 = length(S5);
M10 = load('M_10.mat');      S10 = M10.J_hist;   L10 = length(S10);
M20 = load('M_20.mat');      S20 = M20.J_hist;   L20 = length(S20);
M25 = load('M_25.mat');      S25 = M25.J_hist;   L25 = length(S25);
M50 = load('M_50.mat');      S50 = M50.J_hist;   L50 = length(S50);
M100 = load('M_100.mat');    S100 = M100.J_hist; L100 = length(S100);
M250 = load('M_250.mat');    S250 = M250.J_hist; L250 = length(S250);
M500 = load('M_500.mat');    S500 = M500.J_hist; L500 = length(S500);
t1 = load('T_1.mat');        te1 = t1.telapsed;
v2 = load('T_2.mat');        te2 = v2.telapsed;
v5 = load('T_5.mat');        te5 = v5.telapsed;
v10 = load('T_10.mat');      te10 = v10.telapsed;
v20 = load('T_20.mat');      te20 = v20.telapsed;
v25 = load('T_25.mat');      te25 = v25.telapsed;
v50 = load('T_50.mat');      te50 = v50.telapsed;
v100 = load('T_100.mat');    te100 = v100.telapsed;
v250 = load('T_250.mat');    te250 = v250.telapsed;
v500 = load('T_500.mat');    te500 = v500.telapsed;

v1 = load('V_1.mat');        ve1 = v1.Max_vio;
v2 = load('V_2.mat');        ve2 = v2.Max_vio;
v5 = load('V_5.mat');        ve5 = v5.Max_vio;
v10 = load('V_10.mat');      ve10 = v10.Max_vio;
v20 = load('V_20.mat');      ve20 = v20.Max_vio;
v25 = load('V_25.mat');      ve25 = v25.Max_vio;
v50 = load('V_50.mat');      ve50 = v50.Max_vio;
v100 = load('V_100.mat');    ve100 = v100.Max_vio;
v250 = load('V_250.mat');    ve250 = v250.Max_vio;
v500 = load('V_500.mat');    ve500 = v500.Max_vio;

figure(111);
semilogy(1:5:L1,S1(1:5:L1),'r-o','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L2,S2(1:5:L2),'b-o','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L5,S5(1:5:L5),'g-o','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L10,S10(1:5:L10),'b-x','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L20,S20(1:5:L20),'r-x','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L25,S25(1:5:L25),'g-x','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L50,S50(1:5:L50),'m-x','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L100,S100(1:5:L100),'r-d','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L250,S250(1:5:L250),'g-d','LineWidth',1.0,'MarkerSize',8.0); hold on;
semilogy(1:5:L500,S500(1:5:L500),'c-d','LineWidth',1.0,'MarkerSize',8.0); hold on;
% ylim([1.5 3.5]);
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

%%
figure(333);
ve = [ve1 ve2 ve5 ve10 ve20 ve25 ve50 ve100 ve250 ve500]' * 1e7;
%%%%% final cost---Iter----Computation Time---Max Violation %%%%%
tick = {'Phase: 1','Phase: 2','Phase: 5','Phase: 10','Phase: 20','Phase: 25','Phase: 50','Phase: 100',...
            'Phase: 250','Phase: 500'};
y = [S1(end)   L1   te1  ;
     S2(end)   L2   te2  ;
     S5(end)   L5   te5   ;
     S10(end)  L10  te10  ;
     S20(end)  L20  te20  ;
     S25(end)  L25  te25  ;
     S50(end)  L50  te50  ;
     S100(end) L100 te100 ;
     S250(end) L250 te250 ;
     S500(end) L500 te500 ];
y = [y ve];
b = bar(y,1.0);
ylim([0 65]);
set(gca,'XTickLabel',tick);
legend('Final Cost','Iterations','Computation Time','Constraint Violation ($\times 10^7$)','Interpreter','latex','FontSize',12);
