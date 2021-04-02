clc;clear;
date = '01-Apr-2021';

M1 = load(strcat('M_1_',date,'.mat'));         S1 = M1.J_hist;     J1 = length(S1);   
M2 = load(strcat('M_2_',date,'.mat'));         S2 = M2.J_hist;     J2 = length(S2);
M5 = load(strcat('M_5_',date,'.mat'));         S5 = M5.J_hist;     J5 = length(S5);
M10 = load(strcat('M_10_',date,'.mat'));       S10 = M10.J_hist;     J10 = length(S10);
M20 = load(strcat('M_20_',date,'.mat'));       S20 = M20.J_hist;     J20 = length(S20);
M25 = load(strcat('M_25_',date,'.mat'));       S25 = M25.J_hist;     J25 = length(S25);
M50 = load(strcat('M_50_',date,'.mat'));       S50 = M50.J_hist;     J50 = length(S50);
M100 = load(strcat('M_100_',date,'.mat'));     S100 = M100.J_hist;     J100 = length(S100);
M250 = load(strcat('M_250_',date,'.mat'));     S250 = M250.J_hist;     J250 = length(S250);
M500 = load(strcat('M_500_',date,'.mat'));     S500 = M500.J_hist;     J500 = length(S500);
FJ1 = S1(end);  
FJ2 = S2(end);  
FJ5 = S5(end);  
FJ10 = S10(end); 
FJ20 = S20(end);  
FJ25 = S25(end); 
FJ50 = S50(end);  
FJ100 = S100(end); 
FJ250 = S250(end);  
FJ500 = S500(end); 
% t1 = load(strcat('T_1_',date,'.mat'));        te1 = t1.telapsed;
% t2 = load(strcat('T_2_',date,'.mat'));        te2 = t2.telapsed;
% t5 = load(strcat('T_5_',date,'.mat'));        te5 = t5.telapsed;
% t10 = load(strcat('T_10_',date,'.mat'));      te10 = t10.telapsed;
% t20 = load(strcat('T_20_',date,'.mat'));      te20 = t20.telapsed;
% t25 = load(strcat('T_25_',date,'.mat'));      te25 = t25.telapsed;
% t50 = load(strcat('T_50_',date,'.mat'));      te50 = t50.telapsed;
% t100 = load(strcat('T_100_',date,'.mat'));    te100 = t100.telapsed;
% t250 = load(strcat('T_250_',date,'.mat'));    te250 = t250.telapsed;
% t500 = load(strcat('T_500_',date,'.mat'));    te500 = t500.telapsed;

R_hist1 = load(strcat('R_1_',date,'.mat')); R1 = R_hist1.R_hist;
R_hist2 = load(strcat('R_2_',date,'.mat')); R2 = R_hist2.R_hist;
R_hist5 = load(strcat('R_5_',date,'.mat')); R5 = R_hist5.R_hist;
R_hist10 = load(strcat('R_10_',date,'.mat')); R10 = R_hist10.R_hist;
R_hist20 = load(strcat('R_20_',date,'.mat')); R20 = R_hist20.R_hist;
R_hist25 = load(strcat('R_25_',date,'.mat')); R25 = R_hist25.R_hist;
R_hist50 = load(strcat('R_50_',date,'.mat')); R50 = R_hist50.R_hist;
R_hist100 = load(strcat('R_100_',date,'.mat')); R100 = R_hist100.R_hist;
R_hist250 = load(strcat('R_250_',date,'.mat')); R250 = R_hist250.R_hist;
R_hist500 = load(strcat('R_500_',date,'.mat')); R500 = R_hist500.R_hist;
N = length(S1);

figure(111);
plot(S1,'Color',[0 0.4470 0.7410],'LineStyle','--','LineWidth',2.0); hold on;
plot(S2,'Color',[0.8500 0.3250 0.0980],'LineStyle','-','LineWidth',2.0); hold on;
plot(S5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-','LineWidth',2.0); hold on;
plot(S10,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',2.0); hold on;
plot(S20,'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',2.0); hold on;
plot(S25,'Color',[0.3010 0.7450 0.9330],'LineStyle','-','LineWidth',2.0); hold on;
plot(S50,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',2.0); hold on;
plot(S100,'r','LineWidth',2.0); hold on;
plot(S250,'g','LineWidth',2.0); hold on;
plot(S500,'b','LineWidth',2.0); hold on;
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
Iters = [J1 J2 J5 J10 J20 J25 J50 J100 J250 J500];
Jend =  [FJ1 FJ2 FJ5 FJ10 FJ20 FJ25 FJ50 FJ100 FJ250 FJ500];

set(fig,'defaultAxesColorOrder',[[1 0 0]; [0 0 0]]);
yyaxis left
semilogx(Ms,Iters,'rd--','LineWidth',0.5,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8.0);hold on;
ylabel('Iter','Interpreter','latex','FontSize',15);
yyaxis right
semilogx(Ms,Jend, 'ks--','LineWidth',0.5,'MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8.0);
ylabel('Final Cost','Interpreter','latex','FontSize',15);
title('Different  Shooting Phase','Interpreter','latex','FontSize',15);
xlabel('Shooting Phase','Interpreter','latex','FontSize',15);
grid on;

figure(444);
plot(R1,'Color',[0 0.4470 0.7410],'LineStyle','--','LineWidth',2.0); hold on;
plot(R2,'Color',[0.8500 0.3250 0.0980],'LineStyle','-','LineWidth',2.0); hold on;
plot(R5,'Color',[0.9290 0.6940 0.1250],'LineStyle','-','LineWidth',2.0); hold on;
plot(R10,'Color',[0.4940 0.1840 0.5560],'LineStyle','-','LineWidth',2.0); hold on;
plot(R20,'Color',[0.4660 0.6740 0.1880],'LineStyle','-','LineWidth',2.0); hold on;
plot(R25,'Color',[0.3010 0.7450 0.9330],'LineStyle','-','LineWidth',2.0); hold on;
plot(R50,'Color',[0.6350 0.0780 0.1840],'LineStyle','-','LineWidth',2.0); hold on;
plot(R100,'r','LineWidth',2.0); hold on;
plot(R250,'g','LineWidth',2.0); hold on;
plot(R500,'b','LineWidth',2.0); hold on;
lgd = legend('Shooting Phase: 1','Shooting Phase: 2','Shooting Phase: 5','Shooting Phase: 10','Shooting Phase: 20','Shooting Phase: 25','Shooting Phase: 50','Shooting Phase: 100',...
            'Shooting Phase: 250','Shooting Phase: 500','Interpreter','latex','FontSize',12);
lgd.NumColumns = 2;
ha=gca;
set(ha,'yscale','log');
set(ha,'xscale','log');