%%% single shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
params.dt    = .01;
params.T     = 5;
params.N     = params.T / params.dt;
params.x0    = [0.0; 0.0; 0.0; 0.0];
params.xf    = [3.14; 0.0; 0.0; 0.0];
params.nx    = numel(params.x0);
params.nu    = 1;
params.Q     = diag([0.1 0.1 0.1 0.1])*10;
params.R     =  0.1;
params.Qf    =  diag([0.1 0.1 0.1 0.1])*1200;
params.Rf    = eye(params.nu);
params.Reg_Type = 2;  % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 20.0;
params.umin  = -20.0;
params.Debug = 1;     % 1 -> show details
params.plot = 1;      % 1 -> show plots during optimization
params.Max_iter = 500;
params.stop = 1e-9;
taxis = linspace(0,params.T,params.N);
params.tax = taxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acrobot = Acrobot();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, acrobot, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = ddp_solver(params);
tic
[xbar, ubar, K, du]=solver.Solve(acrobot, cost, params);
toc
% plot
solver.solver_Callback(xbar,ubar,params);
figure(333);hold on;
plot(solver.Jstore,'b-.','LineWidth',2.0);
grid on;

figure(3333);
plot(taxis,ubar(1,:),'r','LineWidth',2.0);
grid on;

figure(444);
plot(taxis,squeeze(K),'LineWidth',2.0);
grid on;