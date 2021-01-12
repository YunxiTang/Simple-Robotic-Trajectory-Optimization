clear;
clc;
ref = load('Ref_Sol.mat');
refsol = ref.sol;

%%
pendubot = Pendubot();
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
FullInput = @(t)([uFit(t)]);
FullState = @(t)([q1Fit(t) q2Fit(t) qd1Fit(t) qd2Fit(t)]);

[k1Fit, ~] = createFit(tn(1:M), K(:,1)');
[k2Fit, ~] = createFit(tn(1:M), K(:,2)');
[k3Fit, ~] = createFit(tn(1:M), K(:,3)');
[k4Fit, ~] = createFit(tn(1:M), K(:,4)');
FullGain = @(t)([k1Fit(t) k2Fit(t) k3Fit(t) k4Fit(t)]);
%%
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

figure(4444);hold on;
plot(q1,q2,'k-.','LineWidth',2.0);

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
%%
x_init = [0.0;0.0;0.0;0.0] + 0.1 * randn(4,1);

Q_ddp = diag([100 100 100 100]);
R_ddp = [0.5];
u_max = 15;
u_min = -15;
Nx = numel(x_init);
Nu = numel(u_max);
num_iter = 2;
alpha = 1.0;
stop_criterion = 1e-3;

% define cost model
cost = QuadraticCostwBarrier(Q_ddp, Q_ddp, R_ddp, u_max, u_min);

%%
N = 10000;
t_scale = linspace(tn(1),tn(end),N);
U_REF = reshape(FullInput(t_scale),[numel(t_scale),Nu]);
X_REF = reshape(FullState(t_scale),[numel(t_scale),Nx]);
K_REF = reshape(FullGain(t_scale),[numel(t_scale),Nx]);

dt = t_scale(2)-t_scale(1);
x_cur = x_init;
x_track = zeros(N,Nx);
U_track = zeros(N,Nu);
CPU_Time = zeros(N,1);
T_upda = 50;
n_h = 80;

%%
% figure;hold on;
tic
for i=1:1:N
    x_track(i,:) = x_cur;
    t_now = t_scale(i);
    if(i > (N - n_h + 1))   % Correction for last n steps
            n_h = n_h - 1;     
    end
    if mod(i,T_upda) == 1
        t_H = t_now + (0:1:n_h-1)*dt;
        x_H = reshape(FullState(t_H),[numel(t_H),Nx]);    % Dim: [nt, Nx] / by rows
        K_H = reshape(FullGain(t_H),[numel(t_H),Nx]);     % Dim: [nt, Nx] / by rows
        u_H = reshape(FullInput(t_H),[numel(t_H),Nu]);    % Dim: [nt, Nu] / by rows
        TSTART = tic;
        [u_rep,K_f,x_rep] = ddp_tracking(x_cur,x_H,u_H,n_h,pendubot,cost,K_H,...
                                         u_max,dt,num_iter,alpha,stop_criterion);
        telapsed = toc(TSTART);
        CPU_Time(i) = telapsed;
        for j=1:1:numel(t_H)
            xrep(j,:) = x_rep{j}';
        end
        for j=1:1:n_h-1
            urep(j,:) = u_rep{j}';
            Kf(j,:) = K_f{j}';
        end
        U_REF(i+1:i+n_h-1,:) = urep(1:n_h-1,:);
        X_REF(i:i+n_h-1,:) = xrep(1:n_h,:);
        K_REF(i+1:i+n_h-1,:) = Kf(1:n_h-1,:);
%         plot(t_H,xrep(:,1),'Color','#4DBEEE','LineWidth',0.2);hold on;
%         plot(t_H,x_H(:,1),'k','LineWidth',2.0);hold on;
    end
%     u_ctrl = u_H(1,:) - K_H(1,:)*(x_cur - x_H(1,:).');
    error = x_cur - X_REF(i,:).';
    u_ctrl = U_REF(i,:) - 0.5 * reshape(K_REF(i,:),[1 4])*error;
    
    U_track(i,:) = u_ctrl;
    x_next = pendubot.rk45(x_cur,u_ctrl,dt);
    if 3 < t_now && t_now <= 3.0 + dt 
        x_next = x_next + [0.0;0.0;0.0;0.0];
    end
    x_cur = x_next;
%     scatter(t_now,x_cur(1),'r.');
%     drawnow;
end
% hold off;
toc
mu = mean(CPU_Time);
sig2 = var(CPU_Time);
%%
figure(2222);
for kk=1:Nx
    subplot(Nx,1,kk);
    plot(t_scale,x_track(:,kk),'LineWidth',2.0);
    title(strcat('$X_{',num2str(kk),'}$'),'Interpreter','latex','FontSize',20);
    grid on;
end
figure(3333);
plot(t_scale,U_track,'LineWidth',2.0);
xlabel('Time [s]','Interpreter','latex','FontSize',15);
title('Control Input','Interpreter','latex','FontSize',20);

figure(4444);
plot(x_track(:,1),x_track(:,2),'LineWidth',2.0);
grid on;
axis equal;
legend('$X_{ref}$','$X_{track}$','Interpreter','latex','FontSize',12);
xlabel('$X_{1}$','Interpreter','latex','FontSize',15);
ylabel('$X_{2}$','Interpreter','latex','FontSize',15);
title('Phase Plot','Interpreter','latex','FontSize',20);

figure(5555);
bar(CPU_Time,25);hold on;
grid on;
title('CPU TIME (NH=70)','Interpreter','latex','FontSize',20);
xlabel('Iteration [-]','Interpreter','latex','FontSize',15);
ylabel('Time [s]','Interpreter','latex','FontSize',15);
ylim([0,T_upda*dt]);