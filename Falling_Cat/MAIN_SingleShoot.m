%%% single shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
params.dt    = .02;
params.T     =  10;
params.N     = params.T / params.dt;
params.x0    = [0.5;0.5;0.5;0.5;0.5; 0.0;0.0;0.0;0.0;0.0];
params.xf    = zeros(10,1);
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([1;1;1;0.1;0.1;1;1;1;0.1;0.1])*10;
params.R     = eye(params.nu)*0.1;
params.Qf    = diag([1;1;1;0.1;0.1;1;1;1;0.1;0.1])*500;
params.Rf    = eye(params.nu);
params.Reg_Type = 2;  % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 10.0;
params.umin  = -10.0;
params.Debug = 1;     % 1 -> show details
params.plot = 1;      % 1 -> show plots during optimization
params.Max_iter = 1000;
params.stop = 1e-9;
taxis = linspace(0,params.T,params.N);
params.tax = taxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fallcat = falling_cat();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, fallcat, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = ddp_solver(params);
tic
[xbar, ubar, K, du]=solver.Solve(fallcat, cost, params);
toc
% plot
solver.solver_Callback(xbar,ubar,params);
hold off;

figure(333);
plot(solver.Jstore,'b-.','LineWidth',2.0);
hold off;
grid on;

figure(3333);
plot(taxis(1:end-1),ubar(1,1:end-1),'r','LineWidth',2.0);hold on;
plot(taxis(1:end-1),ubar(2,1:end-1),'b','LineWidth',2.0);hold off;
legend("$u_L$","$u_R$",'Interpreter','latex','FontSize',12); 
grid on;

%%%%%%%%% animation %%%%%%%%%
figure();
plot(taxis,rad2deg(xbar(1:3,:)),'LineWidth',2.0);