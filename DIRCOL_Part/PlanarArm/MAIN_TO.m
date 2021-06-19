% MAIN.m
% Solve planar arm TO problem
addpath('..\lib\OptimTraj')
clear; clc;

NJ = 5;
model = planar(NJ);
model.gravity = [0 -9.81 0];
rbt = planarArm(model, NJ);

%% params for optimization
params.T = 5.0;
params.dt = 0.01;
params.ulb = -20 * ones(NJ, 1);
params.uub =  20 * ones(NJ, 1);
params.N   = 50;
params.Q   = diag(ones(2*NJ, 1));
params.R   = diag(ones(NJ, 1)) * 0.1;
params.Qf  = diag(ones(2*NJ, 1)) * 10;
params.x0  = [deg2rad(-90);zeros(2*NJ-1, 1)];
params.xf  = [deg2rad(90);zeros(2*NJ-1, 1)];
params.NJ  = NJ;

%% set up function handles
OptDynamics = @(t,x,u)(DynamicsWrapper(t,x,u,@(t,x,u)rbt.Dynamics(t,x,u)));
problem.func.dynamics = OptDynamics;
problem.func.pathObj  = @(t,x,u) pathObjFunc(t,x,u,params);
problem.func.bndObj   = @(t0,x0,tf,xf) bndObjFunc(t0,x0,tf,xf, params);
problem.func.pathCst  = @(t,x,u) pathCstFunc(t, x, u, params);

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
problem.bounds.state.low = 1e3 * -1 * ones(2*NJ, 1);
problem.bounds.state.upp = 1e3 *  1 * ones(2*NJ, 1);

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
t1 = linspace(soln.grid.time(1), soln.grid.time(end), params.T / params.dt + 1);
t2 = linspace(soln.grid.time(1), soln.grid.time(end), params.T / params.dt);
z = soln.interp.state(t1);
u = soln.interp.control(t2);
%% plot solution
figure(1);
plot(t1,rad2deg(z(1:NJ,:)'), 'LineWidth',2.0);
grid on;
title('State');

figure(2);
plot(t2,u', 'LineWidth',2.0);
title('Input');
grid on;
%%
showmotion(model, t1, z(1:NJ,:));
%%
save('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\PlanarArm\PlanarArm_ms_constrained\x_dircol','z');
save('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\PlanarArm\PlanarArm_ms_constrained\u_dircol','u');
save('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\PlanarArm\PlanarArm_ms_constrained\t_dircol','t1');