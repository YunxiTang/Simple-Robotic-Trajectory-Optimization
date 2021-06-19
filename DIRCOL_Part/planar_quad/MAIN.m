% MAIN.m
% Solve planar quadrotor obstacle-avoiding problem
addpath('..\lib\OptimTraj')
clear; clc;
%% load the rbt
rbt = planar_quadrotor();

%% params for optimization
params.T = 6.0;
params.ulb = [ 0.1; 0.1];
params.uub = [ 5.0; 5.0];
params.N = 100;
params.Q                = diag([1 1 1 1 1 1]);
params.R                = diag([0.1 0.1]);
params.Qf               = diag([5 5 5 5 5 5])*10;
params.x0               = [5.0; 2.5; 0.2; 0.0; 0.0; 0.0];
params.xf               = [1.0; 1.5; 0.0; 0.0; 0.0; 0.0];

% obstacles
params.Obstacles = [2.0 3.0 4.0;
                    2.0 1.0 2.0;
                    0.6 0.5 0.6];
%% set up function handles
OptDynamics = @(t,x,u)(DynamicsWrapper(t,x,u,@(t,x,u)rbt.Dynamics(t,x,u)));

problem.func.dynamics = OptDynamics;
problem.func.pathObj = @(t,x,u) pathObjFunc(t,x,u,params);
problem.func.bndObj = @(t0,x0,tf,xf) bndObjFunc(t0,x0,tf,xf, params);
problem.func.pathCst = @(t,x,u) pathCstFunc(t, x, u, params);


% fixed time stepping
problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = params.T;
problem.bounds.finalTime.upp = params.T;

% fixed initial state and final state
problem.bounds.initialState.low = params.x0;
problem.bounds.initialState.upp = params.x0;
problem.bounds.finalState.low = params.xf;
problem.bounds.finalState.upp = params.xf;

% state limitation
problem.bounds.state.low = [-inf; 0;  -pi/6;-inf;-inf;-inf];
problem.bounds.state.upp = [ inf; inf; pi/6; inf; inf; inf];

% control limitation
problem.bounds.control.low = params.ulb;
problem.bounds.control.upp = params.uub;

%% Initial guess at trajectory                          
problem.guess.time = [0,params.T];
problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = [ones(rbt.Nu,1),ones(rbt.Nu,1)];


problem.options.nlpOpt = optimset('Display','iter',...
                                  'MaxFunEvals',2e5,...
                                  'TolFun',1e-7,...
                                  'MaxIter',500,...
                                  'Algorithm','sqp');

problem.options.method = 'trapezoid';
problem.options.trapezoid.nGrid = params.N;
tic
soln = optimTraj(problem);
toc
soln.info

%% Unpack the simulation
t = linspace(soln.grid.time(1), soln.grid.time(end), 300);
z = soln.interp.state(t);
u = soln.interp.control(t);
J = compute_cost(t,z,u,params)
%% plot solution
figure(1);
plot(t,z', 'LineWidth',2.0);
grid on;
title('State');

figure(2);
plot(t,u', 'LineWidth',2.0);
title('Input');
grid on;

%% plot trace
figure(2000);
plot_obstacle(params.Obstacles, 2000);hold on;
plot(params.x0(1), params.x0(2), 'kp', 'MarkerFaceColor', 'b', 'MarkerSize', 15); hold on;
plot(params.xf(1), params.xf(2), 'rh', 'MarkerFaceColor', 'r', 'MarkerSize', 15); hold on;

plot(z(1,:),z(2,:),'r-.','LineWidth',2.0);