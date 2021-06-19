cons_al = load('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\cons_al.mat');
al_const = cons_al.Cons_Vio;
total_const = [al_const solver.Cons_Vio];
figure();
plot(total_const,'k-o','LineWidth',2.0,'MarkerSize',3);hold on;
plot(al_const,'r-o','LineWidth',2.5,'MarkerSize',3.5);
grid on;
title('CartPole Swingup','Interpreter','latex','FontSize',20);
xlabel('Iter.','Interpreter','latex','FontSize',15);
ylabel('Cons. Vio.','Interpreter','latex','FontSize',15);
legend('RLB Stage','AL Stage','Interpreter','latex','FontSize',15);
ha=gca;
set(ha,'xscale','log');