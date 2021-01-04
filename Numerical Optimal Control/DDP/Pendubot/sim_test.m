%%% test simulation
%% post process the reference trajectory and control sequence, 
%% fit to smooth curve or polynomial
clear;
clc;
ref = load('Ref_Sol.mat');
refsol = ref.sol;

Quu = refsol.Q_uu;
Qux = refsol.Q_ux;
Qu = refsol.Q_u;
tn = refsol.t;
N = length(tn);
q1n = zeros(1, N);
q2n = zeros(1, N);
qd1n = zeros(1, N);
qd2n = zeros(1, N);
un = zeros(1, N);
M = N-1;
K = zeros(M,4);

Q_uu = zeros(M,1);
Q_u = zeros(M,1);
Q_ux = zeros(M,4);

for k = 1:N
    q1n(k) = refsol.x{k}(1);
    q2n(k) = refsol.x{k}(2);
    qd1n(k) = refsol.x{k}(3);
    qd2n(k) = refsol.x{k}(4);
    un(k) = refsol.u{k};
end

for i = 1:M
    Q_uu(i,:) = Quu{i}(1);
    Q_u(i,:) = Qu{i}(1);
    Q_ux(i,:) = Qux{i};
    k = -inv(Quu{i}) * Qux{i};
    K(i,:) = -k;
end

[q1Fit, ~] = createFit(tn, q1n);
[q2Fit, ~] = createFit(tn, q2n);
[qd1Fit, ~] = createFit(tn, qd1n);
[qd2Fit, ~] = createFit(tn, qd2n);
[uFit, ~] = createFit(tn, un);

[k1Fit, ~] = createFit(tn(1:M), K(:,1)');
[k2Fit, ~] = createFit(tn(1:M), K(:,2)');
[k3Fit, ~] = createFit(tn(1:M), K(:,3)');
[k4Fit, ~] = createFit(tn(1:M), K(:,4)');


nGrid = 500;
t = linspace(tn(1),tn(end),nGrid);
q1 = q1Fit(t);k1 = k1Fit(t);
q2 = q2Fit(t);k2 = k2Fit(t);
qd1 = qd1Fit(t);k3 = k3Fit(t);
qd2 = qd2Fit(t);k4 = k4Fit(t);
u = uFit(t);

figure(1);
subplot(4,1,1);
plot(t,q1,'r','LineWidth',2.0);hold on;
plot(tn,q1n,'b-*','MarkerSize',3,'LineWidth',1.0);
title('Fitted $q_{1}$','Interpreter','latex','FontSize',15);
subplot(4,1,2);
plot(t,q2,'r','LineWidth',2.0);hold on;
plot(tn,q2n,'b-*','MarkerSize',3,'LineWidth',1.0);
title('Fitted $q_{2}$','Interpreter','latex','FontSize',15);
subplot(4,1,3);
plot(t,qd1,'r','LineWidth',2.0);hold on;
plot(tn,qd1n,'b-*','MarkerSize',3,'LineWidth',1.0);
title('Fitted $\dot{q_{1}}$','Interpreter','latex','FontSize',15);
subplot(4,1,4);
plot(t,qd2,'r','LineWidth',2.0);hold on;
plot(tn,qd2n,'b-*','MarkerSize',3,'LineWidth',1.0);
title('Fitted $\dot{q_{2}}$','Interpreter','latex','FontSize',15);
figure(2);
plot(t,u,'r','LineWidth',2.0);hold on;
plot(tn,un,'bo-','MarkerSize',3,'LineWidth',1.0);
title('Fitted $u$','Interpreter','latex','FontSize',15);


figure(3);
plot(t,k1,t,k2,t,k3,t,k4,'LineWidth',2.0);
legend('$K_1$','$K_2$','$K_3$','$K_4$','Interpreter','latex','FontSize',15);
title('Fitted Gains');

figure(4);
scalex = 1:M;
scaley = [1,2,3,4];
im = imagesc(scalex,scaley,K');
colorbar;
% grid on;
im.AlphaData = .7;
title('Fitted Gains');
%% do a bunch of simulation and plot the resulted trajecory
robot = Pendubot();
testPer = 0.0;
%%
x0 = [0.0;0.0;0.0;0.0] + testPer*randn(4,1);
xf = [pi;0.0;0.0;0.0];
tSpan =[0,tn(end)];
dt = 0.0001;
t_sim = linspace(tSpan(1),tSpan(2),N);
u = @(t,x)(uFit(t)-[k1Fit(t),k2Fit(t),k3Fit(t),k4Fit(t)]*(x-[q1Fit(t);q2Fit(t);qd1Fit(t);qd2Fit(t)]));
dyn = @(t,x)(robot.Dynamics(t,x,u(t,x)));
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
sol = ode45(@(t,x)dyn(t,x),t_sim,x0,opts);

tr = sol.x;
z = sol.y;

ur = [];
for i=1:length(tr)
    ti = t(i);
    zi = z(:,i);
    ur = [ur u(ti,zi)];
end

figure(222);
subplot(4,1,1);hold on;
plot(tr,z(1,:),'LineWidth',1.0);
title('$q_1$','Interpreter','latex','FontSize',15);
subplot(4,1,2);hold on;
plot(tr,z(2,:),'LineWidth',1.0);
title('$q_2$','Interpreter','latex','FontSize',15);

subplot(4,1,3);hold on;
plot(tr,z(3,:),'LineWidth',1.0);
title('$\dot{q}_1$','Interpreter','latex','FontSize',15);
subplot(4,1,4);hold on;
plot(tr,z(4,:),'LineWidth',1.0);
title('$\dot{q}_2$','Interpreter','latex','FontSize',15);
xlabel('Time [s]','Interpreter','latex','FontSize',17);


figure(333);hold on;
plot(tr,ur,'LineWidth',2.0);
title('$u$','Interpreter','latex','FontSize',15);
xlabel('Time [s]','Interpreter','latex','FontSize',17);
%% do a simple animation
[p1,p2] = kinematics(robot.l1, robot.l2, z);
runAnimation(tr, p1, p2);



