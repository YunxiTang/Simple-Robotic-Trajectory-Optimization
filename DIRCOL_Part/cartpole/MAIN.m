% MAIN.m
% Solve the cart-pole swing-up problem

clear; clc;
%% load the cartpole
rbt = CartPole();
Xref = load('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\x_dircol.mat');
Uref = load('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\u_dircol.mat');
Tref = load('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\t_dircol.mat');
xref = Xref.z;
uref = Uref.u;
tref = Tref.t1;
%% params for optimization
params.T = 3.0;
params.dt = 0.01;
params.maxForce = 25.0;
params.N = 10;
params.Q = eye(rbt.Nx);
params.Qf = 5 * eye(rbt.Nx);
params.R = 0.1 * eye(rbt.Nu);
params.x0 = zeros(rbt.Nx, 1);
params.xf = [0.5;pi;0.0;0.0];

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
problem.guess.time = linspace(0, params.T, params.N);
problem.guess.state = [linspace(params.x0(1), params.xf(1), params.N);
                       linspace(params.x0(2), params.xf(2), params.N);
                       linspace(params.x0(3), params.xf(3), params.N);
                       linspace(params.x0(4), params.xf(4), params.N)];
problem.guess.control = zeros(1, params.N);
% problem.guess.time = tref;
% problem.guess.state = xref;
% problem.guess.control = uref;
%%
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
t1 = linspace(soln.grid.time(1), soln.grid.time(end), params.T / params.dt + 1);
t2 = linspace(soln.grid.time(1), soln.grid.time(end), params.T / params.dt);
z = soln.interp.state(t1);
u = soln.interp.control(t2);

%% plot solution
figure(1);
plot(t1,z', 'LineWidth',2.0);
grid on;
title('State');

figure(2);
stairs(t2,u', 'LineWidth',2.0);
title('Input');
grid on;

%%
save('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\x_dircol','z');
save('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\u_dircol','u');
save('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\t_dircol','t1');