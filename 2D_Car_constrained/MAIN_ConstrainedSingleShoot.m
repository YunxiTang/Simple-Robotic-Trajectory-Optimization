%%% Constrained Single Shooting SLQ for Free Robot Dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
log               = 0;       % data log flag  
params.dt         = .01;
params.T          = 10.0;
params.N          = params.T / params.dt;
params.x0         = [1.6;1.5;0.0;0.0];
params.xf         = [2.5;3.5;0.0;0.0];
params.nx         = numel(params.x0);
params.nu         = 2;
params.Q          = diag([0.1 0.1 0.1 0.1])*2;
params.R          = diag([0.2 0.2]);
params.Qf         = diag([500 500 500 500]);
params.Rf         = eye(params.nu);
params.Reg_Type   = 2.0;  % 1->reg of Quu  / 2->reg of Vxx
params.umax       = 5.0;
params.umin       = -5.0;
params.Debug      = 1;     % 1 -> show details
params.plot       = 1;      % 1 -> show plots during optimization
params.Max_iter   = 500;
params.stop       = 1e-9;
taxis             = linspace(0,params.T,params.N);
params.tax        = taxis;
params.qp         = 0;
params.MapNo      = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create Obstacles on the Map %%%
%%%% --------------------------  %%%
%%%% |xc_1 | xc_2 | ... | xc_m|  %%%
%%%% |-----+------+-----+-----|  %%%
%%%% |yc_1 | yc_2 | ... | yc_m|  %%%
%%%% |-----+------+-----+-----|  %%%
%%%% | r1  |  r2  | ... |  rm |  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Obstacles = [0.8 1.5 2.7;
             1.2 2.5 2.5;
             0.5 0.5 0.5];
Constraints = constraint(Obstacles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
car = Car();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[success] = Setup_Functions(params,car,cost,Constraints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = cssddp_solver(Constraints, params);
tic
[xbar, ubar, K, du]=solver.Solve(car, cost, params);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Plotting  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(123);
for kk=1:Constraints.n_ineq
    plot(solver.Lambda(kk,:),'LineWidth',2.0);hold on;
end
grid on;
title('Dual Variables($\lambda$)','Interpreter','latex','FontSize',20);


Primal_residual = zeros(Constraints.n_ineq,params.N);
for m=1:params.N
    Primal_residual(:,m) = Constraints.c_ineq(xbar(:,m), ubar(:,m));
end
figure(234);
for kk=1:Constraints.n_ineq
    plot(Primal_residual(kk,:),'LineWidth',2.0);hold on;
end
plot(zeros(params.N),'r-.','LineWidth',2.0);hold on;
title('Constraint Violation','Interpreter','latex','FontSize',20);
grid on;

figure();
plot(taxis,ubar,'LineWidth',2.0);

