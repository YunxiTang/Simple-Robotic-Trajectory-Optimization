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
params.x0    = [1.0; 1.0; 0.8; 0.0; 0.0; 0.0];
params.xf    = [9.5; 1.0; 0.0; 0.0; 0.0; 0.0];
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([0.1 0.1 0.1 0.1 0.1 0.1])*5;
params.R     = diag([0.1 0.1]);
params.Qf    = diag([50 50 50 50 50 50]);
params.Rf    = eye(params.nu);
params.Reg_Type = 1;  % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 10.0;
params.umin  = -10.0;
params.Debug = 1;     % 1 -> show details
params.plot = 1;      % 1 -> show plots during optimization
params.Max_iter = 500;
params.stop = 1e-9;
taxis = linspace(0,params.T,params.N);
params.tax = taxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
planar_quad = planar_quadrotor();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, planar_quad, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = ddp_solver(params);
tic
[xbar, ubar, K, du]=solver.Solve(planar_quad, cost, params);
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
% figure(789);
% k = 5;
% plot(xbar(1,:),xbar(2,:),'k-.','LineWidth',1.0); hold on;
% planar_quad.animation(taxis,xbar,k,789);
% plot(xbar(1,:),xbar(2,:),'k-.','LineWidth',1.0); hold on;
% plot(xbar(1,1:k:end),xbar(2,1:k:end),'ro','LineWidth',2.0); hold on;