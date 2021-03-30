%%% single shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clear;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
parser.dt    = .02;
parser.T     =  10;
parser.N     = parser.T / parser.dt;
parser.x0    = [1.0; 1.0; 0.8; 0.0; 0.0; 0.0];
parser.xf    = [8.5; 1.0; 0.0; 0.0; 0.0; 0.0];
parser.nx    = numel(parser.x0);
parser.nu    = 2;
parser.Q     = diag([0.1 0.1 0.1 0.1 0.1 0.1])*5;
parser.R     = diag([0.1 0.1]);
parser.Qf    = diag([50 50 50 50 50 50]);
parser.Rf    = eye(parser.nu);
parser.Reg_Type = 1;  % 1->reg of Quu  / 2->reg of Vxx
parser.umax  = 10.0;
parser.umin  = -10.0;
parser.Debug = 1;     % 1 -> show details
parser.plot = 1;      % 1 -> show plots during optimization
parser.Max_iter = 500;
parser.stop = 1e-9;
taxis = linspace(0,parser.T,parser.N);
parser.tax = taxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
planar_quad = planar_quadrotor();
cost = cst_mdl(parser.Q,parser.R,parser.Qf,parser.Rf,parser.umax,parser.umin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(parser, planar_quad, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = ddp_solver(parser);
tic
[xbar, ubar, K, du]=solver.Solve(planar_quad, cost, parser);
toc
% plot
solver.solver_Callback(xbar,ubar,parser);
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