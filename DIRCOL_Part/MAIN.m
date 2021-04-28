% MAIN.m
% Solve the cart-pole swing-up problem

clear; clc;

%% load the cartpole
cartpole = CartPole();

%% params for optimization
params.T = 3.0;
params.maxForce = 35.0;
params.N = 100;
params.Q = eye(cartpole.Nx);
params.R = 0.1 * eye(cartpole.Nu);
params.x0 = zeros(cartpole.Nx, 1);
params.xf = [0.0;pi;0.0;0.0];
%% set up function handles
OptDynamics = @(t,x,u)(DynamicsWrapper(cartpole.Dynamics(t,x,u),x,u));
pathObjFunc = @(t,x,u)((x-params.xf)' * params.Q * (x-params.xf) + u'*params.R*u);

problem.func.dynamics = OptDynamics;
problem.func.pathObj = pathObjFunc;


problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = params.T;
problem.bounds.finalTime.upp = params.T;

problem.bounds.initialState.low = params.x0;
problem.bounds.initialState.upp = params.x0;
problem.bounds.finalState.low = params.xf;
problem.bounds.finalState.upp = params.xf;

problem.bounds.state.low = [-inf;-2*pi;-inf;-inf];
problem.bounds.state.upp = [ inf; 2*pi; inf; inf];

problem.bounds.control.low = -params.maxForce;
problem.bounds.control.upp =  params.maxForce;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time = [0,params.T];
problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = [0,0];

problem.options.nlpOpt = optimset('Display','iter',...
                                  'MaxFunEvals',2e5,...
                                  'TolFun',1e-3);

problem.options.method = 'trapezoid';

tic
soln = optimTraj(problem);
toc

%%%% Unpack the simulation
t = linspace(soln.grid.time(1), soln.grid.time(end), 150);
z = soln.interp.state(t);
u = soln.interp.control(t);