u_g = load('usol1_groundtruth.mat'); ug = u_g.usol;
x_g = load('xsol1_groundtruth.mat'); xg = x_g.xsol;

u_1 = load('usol1.mat'); u1=u_1.usol;
x_1 = load('xsol1.mat'); x1=x_1.xsol;

u_10 = load('usol10.mat'); u10=u_10.usol;
x_10 = load('xsol10.mat'); x10=x_10.xsol;

u_50 = load('usol50.mat'); u50=u_50.usol;
x_50 = load('xsol50.mat'); x50=x_50.xsol;

figure();
plot(ug,'r-','LineWidth',4.0);hold on;
plot(u1,'g-','LineWidth',3.0);hold on;
plot(u10,'b-','LineWidth',2.0);hold on;
plot(u50,'m-','LineWidth',1.0);

figure();
plot(xg(1,:),xg(2,:),'r-','LineWidth',4.0);hold on;
plot(x1(1,:),x1(2,:),'g-','LineWidth',3.0);hold on;
plot(x10(1,:),x10(2,:),'b-','LineWidth',2.0);hold on;
plot(x50(1,:),x50(2,:),'m-','LineWidth',1.0);