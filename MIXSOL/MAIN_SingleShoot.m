%%% single shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
params.dt    = .02;
params.T     = 5;
params.N     = params.T / params.dt;
params.x0    = [0.0;0.0];
params.xf    = [3.14;0.0];
params.nx    = numel(params.x0);
params.nu    = 1;
params.Q     = eye(params.nx);
params.R     = eye(params.nu);
params.Qf    = eye(params.nx);
params.Rf    = eye(params.nu);
params.Reg_Type = 2;  % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 5;
params.umin  = -5;
params.Debug = 1;     % 1 -> show details
params.plot = 0;      % 1 -> show plots during optimization
params.Max_iter = 500;
params.stop = 1e-9;
taxis = linspace(0,params.T,params.N);
params.tax = taxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pendulum = rbt_mdl();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, pendulum, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = ddp_solver(params);
[x,u]=solver.Solve(pendulum, cost, params);
% plot
solver.solver_Callback(x,u,params);