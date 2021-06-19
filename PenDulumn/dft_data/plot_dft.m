dft_10 = load('dft_10.mat');
dft_25 = load('dft_25.mat');
dft_50 = load('dft_50.mat');
dft_100 = load('dft_100.mat');
dft_250 = load('dft_250.mat');

dft10 = dft_10.dftn;
dft25 = dft_25.dftn;
dft50 = dft_50.dftn;
dft100 = dft_100.dftn;
% dft250 = dft_250.dftn;
subplot(2,1,1)
plot(dft10,'k-.s','LineWidth',1.5);hold on;
plot(dft25 ,'r-.d','LineWidth',1.5);hold on;
plot(dft50 ,'g-.+','LineWidth',1.5);hold on;

plot(dft100 ,'b-.*','LineWidth',1.5);hold on;
% plot(dft250/ 250 ,'m-d','LineWidth',1.5);
ha=gca;
% set(ha,'yscale','log');
set(ha,'xscale','log');
ylabel('$|d|_2 $','Interpreter','latex','FontSize',25);
% xlabel('Iteration','Interpreter','latex','FontSize',15);
lgd = legend('M=10','M=25','M=50','M=100','Interpreter','latex','FontSize',15);
lgd.NumColumns = 2;
grid on;