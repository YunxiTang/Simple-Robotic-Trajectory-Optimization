%%% Constrained Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    Date  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_date = date;
log = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
[tref,xref,uref]        = TrajGen();
params.dt               = tref(2) - tref(1);
params.T                = tref(end);
params.N                = size(xref, 2);
params.shooting_phase   = 1;
params.L                = (params.N - 1) / params.shooting_phase + 1;
params.nx               = numel(xref(:,1));
params.nu               = numel(uref(:,1));
%% post-process reference trajectory
params.xref             = cell(params.shooting_phase, 1);
params.tref             = cell(params.shooting_phase, 1);
for i = 1:params.shooting_phase
    params.xref{i} = zeros(params.nx, params.L);
    params.uref{i} = zeros(params.nu, params.L-1);
    params.tref{i} = zeros(params.L, 1);
    
    if 1 < params.shooting_phase
       params.xref{1} = xref(:,1:params.L);
       params.tref{1} = tref(1:params.L);
       params.uref{1} = uref(:,1:params.L-1);
       for k=2:params.shooting_phase                
           params.xref{k} = xref(:, (k-1)*(params.L-1)+(1:params.L));
           params.tref{k} = tref((k-1)*(params.L-1)+(1:params.L));
           params.uref{k} = uref(:,(k-1)*(params.L-1)+(1:params.L-1));
       end 
    else
        params.tref{1} = tref;
        params.xref{1} = xref;
        params.uref{1} = uref;
    end
end

%%
params.Q                = diag([1 1 .1 .1 .1 .1]);
params.R                = diag([0.1 0.1]);
params.Qf               = diag([5 5 5 5 5 5]);
params.Rf               = eye(params.nu);
params.Reg_Type         = 2;     % 1->reg of Quu  / 2->reg of Vxx
params.umax             = 5.0;
params.umin             = 0.1;
params.Debug            = 1;     % 1 -> show details
params.plot             = 1;     % 1 -> show plots during optimization
params.Max_iter         = 500;
params.stop             = 1e-3;
params.qp               = 0;
params.clamp            = 1;
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
Obstacles = 0*[2.0 3.0 4.0;
             2.0 1.0 2.0;
             0.6 0.5 0.6];

Constraints = constraint(Obstacles, params.dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
planar_quad = planar_quadrotor();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, planar_quad, cost, Constraints);

%% ..................................................
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = cmsddp_solver(Constraints, params);

tstart = tic;
[xsol, usol, Ksol, Lambdasol, dft, xbar, ubar] = solver.Solve(planar_quad,cost,params);
telapsed = toc(tstart)

figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0,'MarkerSize',3);
ha=gca;
set(ha,'yscale','log');
set(ha,'xscale','log');
grid on;
J_hist = solver.Jstore;

%%
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

%%%
figure(2000);
Constraints.plot_obstacle(2000);hold on;
% plot(params.x0(1), params.x0(2), 'kp', 'MarkerFaceColor', 'b', 'MarkerSize', 15); hold on;
% plot(params.xf(1), params.xf(2), 'rh', 'MarkerFaceColor', 'r', 'MarkerSize', 15); hold on;
h = plot(xsol(1,:),xsol(2,:),'r-.','LineWidth',2.0); hold on;
plot(xref(1,:), xref(2,:), 'b--','LineWidth',2.0)

grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% data logging %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if log == 1
    file_name1 = strcat('.\data\T_', ... 
                         num2str(params.shooting_phase));
    file_name2 = strcat('.\data\M_', ...
                         num2str(params.shooting_phase));
    file_name3 = strcat('.\data\V_', ...
                         num2str(params.shooting_phase));
    save(file_name1,'telapsed');
    save(file_name2,'J_hist');
    save(file_name3,'Max_vio');
end

% save coarse solution
filename_x = '.\x_al';
save(filename_x,'xbar');
filename_u = '.\u_al';
save(filename_u,'ubar');