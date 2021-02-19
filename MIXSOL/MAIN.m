%%% mutiple shooting SLQ for free robot dynamics 
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
params.umax  = 5;
params.umin  = -5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);
pendulum = rbt_mdl();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, pendulum, cost);


