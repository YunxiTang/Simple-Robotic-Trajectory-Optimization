%%% Test the obtained trajectory and controller
clear;
clc;
% rng('default');
%% Load optimized solutions
tsol = load('.\Opt_Solution\t_sol.mat');
xsol = load('.\Opt_Solution\x_sol.mat');
usol = load('.\Opt_Solution\u_sol.mat');
Ksol = load('.\Opt_Solution\K_sol.mat');

tref = tsol.t;
xref = xsol.xsol;
uref = usol.usol;
Kref = Ksol.Ksol;

%% Build a robot
% planar_quad = me_planar_quadrotor();
planar_quad = planar_quadrotor();
params.dt = .01;
params.T  = 5.0;
params.x0 = xref(:,1);
params.xf = [1.0; 1.5; 0.0; 0.0; 0.0; 0.0];
params.umax = 5.0;
params.umin = 0.1;
nx = 6;
nu = 2;
%% simulation
dt = 0.0001;
t_hist =  0:dt:params.T;
N = length(t_hist);
x_hist = zeros(nx, N);
u_hist = zeros(nu, N);
x = params.x0 + 0.1 * randn(6,1);
x_hist(:,1) = x;
for i = 1:N-1
    ti = t_hist(i);
    k = max(get_idx(ti,tref),1);
    xi_b = xref(:,k);
    ui_b = uref(:,k);
    Ki_b = Kref(:,:,k);
    u = clamp(ui_b + Ki_b * (x - xi_b), params.umax, params.umin);
    x_next = planar_quad.rk(x, u, dt);
    x = x_next;
    x_hist(:,i+1) = x;
    u_hist(:,i) = u;
end

%% plot
Obstacles = [2.0 3.0 4.0;
             2.0 1.0 2.0;
             0.6 0.5 0.6];
figure(2000);
k = 25;
plot_obstacle(Obstacles-[0 0 0;0 0 0;0.00 0.00 0.00], 2000);
draw_circle(params.x0(1), params.x0(2),0.1/3);
hold on;
% plot(params.x0(1), params.x0(2), 'bp', 'MarkerFaceColor', 'b', 'MarkerSize', 15); hold on;
plot(params.xf(1), params.xf(2), 'rh', 'MarkerFaceColor', 'r', 'MarkerSize', 15); hold on;
plot(x_hist(1,:), x_hist(2,:),'-.','LineWidth',2.0);
hold on;

figure(1);
stairs(t_hist, u_hist(1,:),'LineWidth',2.0);hold on;
stairs(t_hist, u_hist(2,:),'LineWidth',2.0);hold off;

figure(2);
plot(t_hist, x_hist,'LineWidth',2.0);hold on;

% figure(3);
% plot(t_hist, x_hist - xref, 'LineWidth',2.0);hold on;
%%
function k = get_idx(t,tref)
    dt = tref(2) - tref(1);
    n = floor(t/dt);
    res = t - n*dt;
    if res < dt / 2
        k = n;
    else
        k = n + 1;
    end
end

function [] = draw_circle(x,y,r)
    xc = x;
    yc = y;
    d = r;
    ds = 0:0.1:2*pi;
    xx = xc + d * cos(ds);
    yy = yc + d * sin(ds);
    plot(xx, yy, 'm--','LineWidth',0.5);
    hold on;
end

function [u_sat] = clamp(u, u_max, u_min)
    nu = numel(u);
    ub = u_max * ones(nu, 1);
    lb = u_min * ones(nu, 1);
    u_sat = max(lb, min(ub, u));
end