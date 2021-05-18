% MAIN.m
% Solve the cart-pole swing-up problem

clear; clc;
%% load the cartpole
rbt = CartPole();

%% params for optimization
params.T = 3.0;
params.maxForce = 25.0;
params.N = 150;
params.Q = eye(rbt.Nx);
params.Qf = 5 * eye(rbt.Nx);
params.R = 0.1 * eye(rbt.Nu);
params.x0 = zeros(rbt.Nx, 1);
params.xf = [0.0;pi;0.0;0.0];

%% set up function handles
OptDynamics = @(t,x,u)(DynamicsWrapper(t,x,u,@(t,x,u)rbt.Dynamics(t,x,u)));

problem.func.dynamics = OptDynamics;
problem.func.pathObj = @(t,x,u) pathObjFunc(t,x,u,params);
problem.func.bndObj = @(t0,x0,tf,xf) bndObjFunc(t0,x0,tf,xf, params);

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
problem.bounds.state.low = [-0.8;-1.5*pi;-inf;-inf];
problem.bounds.state.upp = [ 0.8; 1.5*pi; inf; inf];

% control limitation
problem.bounds.control.low = -params.maxForce;
problem.bounds.control.upp =  params.maxForce;

%% Initial guess at trajectory                          
problem.guess.time = [0,params.T];
problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = [0,0];


problem.options.nlpOpt = optimset('Display','iter',...
                                  'MaxFunEvals',5e5,...
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

%% plot solution
figure(1);
plot(t,z', 'LineWidth',2.0);
grid on;
title('State');

figure(2);
plot(t,u', 'LineWidth',2.0);
title('Input');
grid on;