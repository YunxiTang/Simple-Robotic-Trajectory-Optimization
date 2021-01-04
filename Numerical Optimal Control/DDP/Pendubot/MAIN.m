%%% Differential Dynamic Programming for Dynamical System, including both
%%% trajectory optimization and feedback controller synthesis
clc;
clear;
% close all;
rng(0);

%% problem definition
% Initial state
x_0 = [0.0; 0.0; 0.0; 0.0];

% Target state
x_star = [3.14; 0.0; 0.0; 0.0];

% Time horizon
t_f = 10.0;

% Number of states along trajectory
N = floor(t_f ./ 0.025);

% Maximum magnitude of control
u_max = [10.0];

% Initialize dynamics
fprintf("Initializing dynamics...\n")

dyn = Pendubot();

% Initialize cost
fprintf("Initializing quadratic cost function...\n")
Q_f = 140 * eye(4);
Q = [1000  0  0     0;
     0     1000  0     0;
     0     0  1000  0;
     0     0     0  1000];
R = 0.5;
cost = QuadraticCost_Barrier(50*Q, Q/60, R, u_max, -u_max);

% Max Number of DDP iterations
num_iter = 500;

% Stop criterion
stop_criterion = 1e-7;

% DDP learning rate
alpha = 0.01;

% Video framerate
fps = 30;

%% Execution of DDP

fprintf("executing DDP...");

tic;
sol = DDP(x_0, x_star, t_f, N, dyn, cost, u_max, num_iter, alpha,stop_criterion);
toc;
save('Ref_Sol','sol');
%% Begin post-processing of solution

if sol.error == 1
    fprintf("DIVERGENCE ERROR: try decreasing learning rate\n");
    return;
end

q1 = zeros(1, length(sol.x));
q2 = zeros(1, length(sol.x));
qd1 = zeros(1, length(sol.x));
qd2 = zeros(1, length(sol.x));
u = zeros(1, length(sol.x));

for k = 1:N
    q1(k) = sol.x{k}(1);
    q2(k) = sol.x{k}(2);
    qd1(k) = sol.x{k}(3);
    qd2(k) = sol.x{k}(4);
    u(k) = sol.u{k};
end
save('Norm_u.mat','u');
%% Plot trajectories
% Plot cart position history
figure(1);hold on;
pbaspect([5 3 1])
hold on;
plot(sol.t, q1, "Linewidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("$q_{1}$ [rad]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot cart velocity history
figure(2);hold on;
pbaspect([5 3 1])
hold on;
plot(sol.t, qd1, "Linewidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("$q_{d1}$ [m/s]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot pole angle history
figure(3);hold on;
pbaspect([5 3 1])
hold on;
plot(sol.t, q2, "Linewidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("$q_{2}$ [rad]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot pole anglular velocity history
figure(4);hold on;
pbaspect([5 3 1])
hold on;
plot(sol.t, qd2, "Linewidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("$q_{d2}$ [rad/s]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot trajectory in partial state space
figure(5);hold on;
pbaspect([5 3 1])
hold on;
plot(q1, q2,'k-.', 'Linewidth', 2);
plot(x_0(1), x_0(2), 'o', "MarkerFaceColor", "blue", ...
                          "MarkerEdgeColor", "blue");
plot(x_star(1), x_star(2), 'o', "MarkerFaceColor", "green", ...
                                "MarkerEdgeColor", "green")
grid on;
xlabel("$q_{1}$ [m]", "Interpreter", "latex", "FontSize", 20);
ylabel("$q_{2}$ [rad]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot control sequence
figure(6);hold on;
pbaspect([5 3 1])
hold on;
plot(sol.t, u, "LineWidth", 2);
grid on;
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
ylabel("Control Input [N]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";

% Plot cost function vs iteration
figure(7);hold on;
pbaspect([5 3 1])
hold on;
plot(1:length(sol.J), sol.J, "LineWidth", 2);
grid on;
xlabel("DDP Iteration [-]", "Interpreter", "latex", "FontSize", 20);
ylabel("Cost Function [-]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";
ax.YScale = 'log';

% Plot control energy vs iteration
figure(8);hold on;
pbaspect([5 3 1])
hold on;
plot(1:length(sol.E), sol.E, "LineWidth", 2);
grid on;
xlabel("DDP Iteration [-]", "Interpreter", "latex", "FontSize", 20);
ylabel("Control Energy Usage [$\rm{N}^{2}\rm{s}$]", "Interpreter", "latex", "FontSize", 20);
ax = gca();
ax.FontSize = 16;
ax.TickLabelInterpreter = "latex";
ax.YScale = 'log';

%% plot gain and some matrix
compute_feedback_gain;

%% run animation

for i=1:length(q1)
    x(1,:) = q1;
    x(2,:) = q2;
    x(3,:) = qd1;
    x(4,:) = qd2;
end
[p1,p2] = kinematics(dyn.l1, dyn.l2, x);
% runAnimation(sol.t, p1, p2);