%%% Differential Dynamic Programming(SLQ/iLQR) w/ MPC for Dynamical System, including both
%%% trajectory generation and trajectory stablization

clc;clear;
param = struct('l1',1.0,'l2',1.0,...
               'm1',1.0,'m2',1.0,...
               'lc1',1.5,'lc2',0.5,...
               'b1',0.1,'b2',0.1,...
               'I1',0.33,'I2',0.33,...
               'g',9.81);

acrobot = Acrobot_mdl(param);

%% load ref trajs
ref = load('Ref_Sol.mat');refsol = ref.sol;
Quu = refsol.Q_uu; Qux = refsol.Q_ux;
Qu = refsol.Q_u;   tn = refsol.t;
N = length(tn);
q1n = zeros(1, N);   q2n = zeros(1, N);
qd1n = zeros(1, N);  qd2n = zeros(1, N);
un = zeros(1, N);
M = N-1;
K = zeros(M,4);

Q_uu = zeros(M,1); Q_u = zeros(M,1); Q_ux = zeros(M,4);
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
%% 
x_ref = [q1n; q2n; qd1n; qd2n];
u_ref = un(:,1:end-1);
Q_lqr = 1500*eye(4);
R_lqr = [100];

Q_ddp = 15000*diag([1 1 1 1]);
R_ddp = [10];

dt = tn(end) / N;

x_init = [0.0;0.5;0.0;0.0];
x_now = x_init;
u_now = un(1);
n_steps = 10;
u_max = 15;
u_min = -15;
num_iter = 50;
alpha = 1.0;
stop_criterion = 1e-4;

control=[];  
T = [0];
x_track = zeros(4,N);
CTR_Step = 15;
% define cost model
cost = QuadraticCostwBarrier(Q_ddp, Q_ddp, R_ddp, u_max, u_min);
%%
for tt=1:1:M
    T = [T tt*dt];
    x_track(:,tt) = x_now;
    
    % when getting close to the last n_steps
    if(tt > (M - n_steps + 1))   
        n_steps = n_steps - 1;     
    end
    x_des_slq = x_ref(:,(tt+1):(tt+n_steps));
    u_des_slq = u_ref(:,(tt):(tt+n_steps-1));
    [A_0, B_0] = acrobot.linearize(x_des_slq(:,1), u_des_slq(:,1));
    [K0_lqr,S0_lqr,~] = lqr(A_0, B_0, Q_lqr, R_lqr);
    x_tf = x_des_slq(:,end);
    u_tf = u_des_slq(:,end);
    [A_f, B_f] = acrobot.linearize(x_tf, u_tf);
    [Kf_lqr,Sf_lqr,~] = lqr(A_f, B_f, Q_lqr, R_lqr);
    cost.update_Qf(Sf_lqr);
    tic
    [u_r,Kk,x_rep] = ddp_solver(x_now, x_des_slq, u_des_slq, n_steps, acrobot, cost, K(tt,:),...
                                u_max, dt, num_iter, alpha, stop_criterion);
    toc
    uu = u_r - Kk*(x_now-x_rep{1});
    x_next = acrobot.rk45(x_now, uu, dt);
    control =[control uu];
    x_now = x_next;
    
    for p=1:length(x_rep)
        xrep(p) = x_rep{p}(1);
    end
    plot(tn(tt+1:tt+n_steps),x_des_slq(1,:),'-.','LineWidth',2.0);hold on;
    plot(tn(tt+1:tt+n_steps),x_des_slq(2,:),'b-.','LineWidth',2.0);hold on;
    plot(tn(tt+1:tt+n_steps),xrep,'LineWidth',4.0);
    scatter(tt*dt,x_now(1),'g.');
    drawnow;
    pause;
end
%%
figure(111);
plot(T,x_track(1,:),'LineWidth',2.0);hold on;
plot(tn,x_ref(1,:),'-.','LineWidth',2.0);
figure(222);
plot(tn(1:399),u_ref,'b-.','LineWidth',4.0);hold on;
plot(tn(1:399),control,'m','LineWidth',2.0);