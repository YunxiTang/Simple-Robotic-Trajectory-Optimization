%%% multiple forward test
clc;
clear;

%% problem definition
% Initial state
x_0 = [0.0; 0.0; 0.0; 0.0];
% Maximum magnitude of control
u_max = [10.0];
Nx = numel(x_0);
Nu = numel(u_max);

% Target state
x_star = [pi; 0.0; 0.0; 0.0];
% Time horizon
t_0 = 0.0;
t_f = 10.0;
% Number of states along trajectory
N = floor(t_f ./ 0.01)+1;

% Initialize dynamics
pendubot = Pendubot();

% Initialize cost
Q_f = 1400 * eye(4);
Q = 140 * eye(4);
R = 5;
cost = QuadraticCostwBarrier(Q_f, Q, R, u_max, -u_max);

% Max Number of DDP iterations
num_iter = 500;

% Stop criterion
stop_criterion = 1e-7;

% DDP learning rate
alpha = 1.0;

%% break into M intervals
S = 25;
M = N-1;
L = M / S;
dt = (t_f - t_0) / M;
%% x guess
x_interval_guess = [linspace(x_0(1),x_star(1),S+1);
                    linspace(x_0(2),x_star(2),S+1);
                    linspace(x_0(3),x_star(3),S+1);
                    linspace(x_0(4),x_star(4),S+1)];
u_guess = zeros(Nu, M);
t_interval = linspace(t_0,t_f,S+1);
%% first time of forward similation (multiple shooting)
k = 1;
X = zeros(S, Nx, L+1);
U = zeros(S, Nu, L);
T = zeros(S, Nu, L+1);
U_guess = U;
t_now = 0.0;

for i=1:S
    xs_0 = x_interval_guess(:,i);
    xc = xs_0;
    u_now = U_guess(i,:,1);
    X(i,:,1) = xc;
    T(i,:,1) = t_now;
    for j=1:L
        t_now = k * dt;
%         u_now = U_guess(i,:,j)+1.0;
        u_now = 0.00*sin(t_now);
        x_next = pendubot.rk45(xc, u_now, dt); 
        X(i,:,j+1) = x_next;
        T(i,:,j+1) = t_now;
        U(i,:,j) = u_now;
        xc = x_next;
        k = k + 1;
    end
end
%%
figure(111);
h1 = plot(t_interval,x_interval_guess(1,:),'bo','LineWidth',0.8);hold on;
for k=1:S
    h2 = plot(squeeze(T(k,1,:)),squeeze(X(k,1,:)),'k','LineWidth',2.0);hold on;
    h3 = plot(squeeze(T(k,1,1:end-1)),squeeze(U(k,1,:)),'r','LineWidth',2.0);hold on;
end
grid on;
legend([h1,h2,h3],'$Initial\;Guess$','$Mutiple\;Shooting$','$Control\;Input$',...
       'Interpreter','latex','Fontsize',12);
title('$x_1$','Interpreter','latex','Fontsize',18);
xlabel('$Time [s]$','Interpreter','latex','Fontsize',15);
hold off;


