clc;clear;

M1 = load(strcat('M_1_RLB','.mat'));         S1 = M1.J_hist;     J1 = length(S1);   
M2 = load(strcat('M_1_RLB_boxQP','.mat'));   S2 = M2.J_hist;     J2 = length(S2);

t1 = load(strcat('T_1_RLB','.mat'));        te1 = t1.telapsed;
t2 = load(strcat('T_1_RLB_boxQP','.mat'));  te2 = t2.telapsed;

figure(111);
plot(S1,'r-','LineWidth',2.0); hold on;
plot(S2,'k-','LineStyle','-','LineWidth',2.0); 
ha=gca;
set(ha,'yscale','log');set(ha,'xscale','log');
grid on;
title('Cost v.s. Iterations','Interpreter','latex','FontSize',15);
xlabel('Iter','Interpreter','latex','FontSize',15);
ylabel('Cost','Interpreter','latex','FontSize',15);
lgd = legend('$RLB$','$RLB+boxQP$','Interpreter','latex','FontSize',12);