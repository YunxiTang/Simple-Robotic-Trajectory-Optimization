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
NJ = 3;
params.dt               =  .01;
params.T                =  5.0;
params.N                = params.T / params.dt;
params.shooting_phase   = 50;
params.x0               = [deg2rad(-90);zeros(2*NJ-1, 1)];
params.xf               = [deg2rad(90);zeros(2*NJ-1, 1)];
params.nx               = numel(params.x0);
params.nu               = NJ;
params.Q                = diag(ones(2*NJ, 1));
params.R                = diag(ones(NJ, 1)) * 0.1;
params.Qf               = diag(ones(2*NJ, 1)) * 10;
params.Rf               = eye(params.nu);
params.Reg_Type         = 2;     % 1->reg of Quu  / 2->reg of Vxx
params.umax             = 25;
params.umin             = -25;
params.Debug            = 1;     % 1 -> show details
params.plot             = 1;     % 1 -> show plots during optimization
params.Max_iter         = 200;
params.stop             = 1e-7;
params.qp               = 0;
params.clamp            = 1;
params.warm_start       = 0;
params.loop             = 0;
nt                      = params.T / params.shooting_phase;
tax                     = cell(params.shooting_phase,1);
for i=1:params.shooting_phase
    tax{i}=linspace((i-1)*nt,i*nt,nt/params.dt+1);
end
params.t = tax;
params.L = params.N / params.shooting_phase + 1; 
params.MapNo = 1;

%% post-process reference trajectory
Xref = load('.\x_dircol.mat');
Uref = load('.\u_dircol.mat');
Tref = load('.\t_dircol.mat');
xref = Xref.z;
uref = Uref.u;
tref = Tref.t1;
params.xref = cell(params.shooting_phase, 1);
params.tref = cell(params.shooting_phase, 1);
for i = 1:params.shooting_phase
    params.xref{i} = zeros(params.nx, params.L);
    params.uref{i} = zeros(params.nu, params.L-1);
    params.tref{i} = zeros(1, params.L);
    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Add Path Constraints     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_const_func = @(x,u)([]);
Npath_ineq = numel(path_const_func(zeros(params.nx),zeros(params.nu)));
path_cons = path_constraint(path_const_func, Npath_ineq, params.dt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Add Final Constraints    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_const_func = @(xf  )([]);
Nfinal_ineq = numel(final_const_func(zeros(params.nx)));
final_cons = final_constraint(final_const_func,Nfinal_ineq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = planar(NJ);
% model.gravity = [0 -9.81 0];
rbt = planararm(model, NJ);
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, rbt, cost, path_cons, final_cons);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = cmsddp_solver(path_cons, final_cons, params);

tstart = tic;
[xsol, usol, Ksol, dft, xbar, ubar] = solver.Solve(rbt,cost,path_const_func,final_const_func,params);
telapsed = toc(tstart)
solver.iter
%%
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
%% save constraint violation
figure(666);
plot(solver.Cons_Vio,'k-o','LineWidth',2.0,'MarkerSize',3);
filename_cons = '.\cons_al';
Cons_Vio = solver.Cons_Vio;
save(filename_cons,'Cons_Vio');

%%
figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0,'MarkerSize',3);
ha=gca;
set(ha,'yscale','log');
set(ha,'xscale','log');
grid on;
J_hist = solver.Jstore;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  plot data   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% figure(123);
% for kk=1:path_cons.n_ineq
%     plot(Lambdasol(kk,:),'LineWidth',2.0);hold on;
% end
% grid on;
% title('Dual Variables($\lambda$)','Interpreter','latex','FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0.0:params.dt:params.T;
figure();
plot(t, xsol,'LineWidth',2.0); 
legend("$x$","$\theta$","$\dot{x}$","$\dot{\theta}$",'Interpreter','latex','FontSize',12);
grid on;
figure();
plot(t(1:end-1), usol,'LineWidth',2.0);
legend("$u_{F}$",'Interpreter','latex','FontSize',12);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

real_cost = compute_cost(xsol,usol,cost,params);

Primal_residual = zeros(path_cons.n_ineq,params.N);
Vio_Norm = zeros(params.N,1);
for m=1:params.N
    h = path_cons.c(xsol(:,m), usol(:,m));
    Primal_residual(:,m) = h;
    Vio_Norm(m) = norm(max(h,0));
end
Max_vio = max(Vio_Norm);
figure(234);
for kk=1:path_cons.n_ineq
    plot(Primal_residual(kk,:),'LineWidth',2.0);hold on;
end
plot(zeros(params.N),'r-.','LineWidth',2.0);hold on;
title('Constraint Violation','Interpreter','latex','FontSize',20);
grid on;

%%
showmotion(model, t, xsol(1:NJ,:));