du_10 = load('du_10.mat');
du_25 = load('du_25.mat');
du_50 = load('du_50.mat');
du_100 = load('du_100.mat');
% du_250 = load('du_250.mat');

du10 = du_10.dun;
du25 = du_25.dun;
du50 = du_50.dun;
du100 = du_100.dun;
% dft250 = dft_250.dftn;
subplot(2,1,2)
plot(du10,'k-.s','LineWidth',1.5);hold on;
plot(du25 ,'r-.d','LineWidth',1.5);hold on;
plot(du50 ,'g-.+','LineWidth',1.5);hold on;

plot(du100 ,'b-.*','LineWidth',1.5);hold on;
% plot(dft250/ 250 ,'m-d','LineWidth',1.5);
ha=gca;
set(ha,'yscale','log');
% set(ha,'xscale','log');
ylabel('$|du|_2 $','Interpreter','latex','FontSize',25);
xlabel('Iteration','Interpreter','latex','FontSize',15);
lgd = legend('M=10','M=25','M=50','M=100','Interpreter','latex','FontSize',15);
lgd.NumColumns = 2;
grid on;