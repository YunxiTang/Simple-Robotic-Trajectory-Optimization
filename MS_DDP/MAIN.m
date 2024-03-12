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
t_f = 8.0;
% Number of states along trajectory
N = floor(t_f ./ 0.01) + 1;

% Initialize dynamics
pendubot = Pendubot();

% Initialize cost
Q_f = 1000 * eye(4);
Q = 4 * eye(4);
R = 1;
cost = QuadraticCostwBarrier(Q_f, Q, R, u_max, -u_max);

% Max Number of DDP iterations
num_iter = 500;

% Stop criterion
stop_criterion = 1e-7;

% DDP learning rate
alpha = 0.1;

% TRAIL: to store iterated trajectory trails
TRAIL = sol_mdl(num_iter);

%% break into M intervals
S = 20;
M = N-1;
L = M / S;
dt = (t_f - t_0) / M;

tic;
sol = ms_ddp(x_0,x_star,pendubot,cost,u_max,N,M,S,dt,num_iter,alpha,stop_criterion,TRAIL);
toc;