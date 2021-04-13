%%% Constrained Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    Date  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_date = date;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
log = 0;
params.dt               =  .02;
params.T                =  10.0;
params.N                = params.T / params.dt;
params.shooting_phase   = 100;
params.x0               = [0.0;0.0;0.0;0.0];
params.xf               = [2.8;2.8;pi/2;0.0];
params.nx               = numel(params.x0);
params.nu               = 2;
params.Q                = diag([0.1 0.1 0.1 0.1]);
params.R                = diag([0.1 0.1]);
params.Qf               = diag([50 50 50 50])*10;
params.Rf               = eye(params.nu);
params.Reg_Type         = 1;  % 1->reg of Quu  / 2->reg of Vxx
params.umax             = 3.5;
params.umin             = -3.5;
params.Debug            = 1;     % 1 -> show details
params.plot             = 1;     % 1 -> show plots during optimization
params.Max_iter         = 500;
params.stop             = 1e-9;
params.qp               = 1;
params.clamp            = 0;
nt                      = params.T / params.shooting_phase;
tax                     = cell(params.shooting_phase,1);
for i=1:params.shooting_phase
    tax{i}=linspace((i-1)*nt,i*nt,nt/params.dt+1);
end
params.t = tax;

params.MapNo = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create Obstacles on the Map %%%
%%%% --------------------------  %%%
%%%% |xc_1 | xc_2 | ... | xc_m|  %%%
%%%% |-----+------+-----+-----|  %%%
%%%% |yc_1 | yc_2 | ... | yc_m|  %%%
%%%% |-----+------+-----+-----|  %%%
%%%% | r1  |  r2  | ... |  rm |  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obstacles = [0.0  1.0 0.9 0.5 1.8  1.9;
%              1.0  3.0 1.5 2.3 2.8  2.0;
%              0.3  0.3 0.3 0.4 0.25 0.2];
% Obstacles = [0.6 0.5 2.0 3.0;
%              0.4 1.0 2.0 2.4;
%              0.1 0.3 0.5 0.2];
Obstacles = [1.5 ;
             1.5 ;
             0.7 ];
Constraints = constraint(Obstacles, params.dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
car = Car();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, car, cost, Constraints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = cmsddp_solver(Constraints, params);

tstart = tic;
[xsol, usol, Ksol, Lambdasol, dft] = solver.Solve(car,cost,params);
telapsed = toc(tstart)

figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
J_hist = solver.Jstore;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  plot data   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(123);
for kk=1:Constraints.n_ineq
    plot(Lambdasol(kk,:),'LineWidth',2.0);hold on;
end
grid on;
title('Dual Variables($\lambda$)','Interpreter','latex','FontSize',20);

Primal_residual = zeros(Constraints.n_ineq,params.N);
Vio_Norm = zeros(params.N,1);
for m=1:params.N
    h = Constraints.c_ineq(xsol(:,m), usol(:,m));
    Primal_residual(:,m) = h;
    Vio_Norm(m) = norm(max(h,0));
end
Max_vio = max(Vio_Norm);
figure(234);
for kk=1:Constraints.n_ineq
    plot(Primal_residual(kk,:),'LineWidth',2.0);hold on;
end
plot(zeros(params.N),'r-.','LineWidth',2.0);hold on;
title('Constraint Violation','Interpreter','latex','FontSize',20);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0.0:params.dt:params.T;
figure();
plot(t, xsol,'LineWidth',2.0); 
legend("$x$","$y$","$\theta$","$v$",'Interpreter','latex','FontSize',12);
grid on;
figure();
plot(t(1:end-1), usol,'LineWidth',2.0);
legend("$u_{\theta}$","$u_v$",'Interpreter','latex','FontSize',12);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% data logging %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if log == 1
    file_name1 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\2D_Car_ms_constrained\data\T_', ... 
                         num2str(params.shooting_phase));
    file_name2 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\2D_Car_ms_constrained\data\M_', ...
                         num2str(params.shooting_phase));
    file_name3 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\2D_Car_ms_constrained\data\V_', ...
                         num2str(params.shooting_phase));
    save(file_name1,'telapsed');
    save(file_name2,'J_hist');
    save(file_name3,'Max_vio');
end