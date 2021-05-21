%%% Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
% run('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\MAIN_CMShoot.m');
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    Date  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_date = 'RLB';
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 
logger = 0;
params.dt               =  0.01;
params.T                =  3.0;
params.N                = params.T / params.dt;
params.shooting_phase   = 5;
params.x0               = [0.0;0.0;0.0;0.0];
params.xf               = [0.0;pi;0.0;0.0];
params.nx               = numel(params.x0);
params.nu               = 1;
params.Q                = diag([1 1 1 1]);
params.R                = diag([0.1]);
params.Qf               = diag([5 5 5 5]);
params.Rf               = eye(params.nu);
params.Reg_Type         = 2;        % 1->reg of Quu  / 2->reg of Vxx
params.umax             = 25;
params.umin             = -25;
params.Debug = 1;           % 1 -> show details
params.plot = 1;            % 1 -> show plots during optimization
params.Max_iter = 20;
params.stop = 1e-6;
params.qp = 0;
params.warm_start = 1;
nt = params.T / params.shooting_phase;
tax = cell(params.shooting_phase,1);
for i=1:params.shooting_phase
    tax{i}=linspace((i-1)*nt,i*nt,nt/params.dt+1);
end
params.t = tax;
t = 0.0:params.dt:params.T;

if params.warm_start == 1
   x_al = load('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\x_al.mat');
   u_al = load('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\CartPole\CartPole_ms_constrained\u_al.mat');
   params.x_warm = x_al.xbar;
   params.u_warm = u_al.ubar;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot,constraints & cost mdl %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% build a robot
rbt = rbt_mdl();
% cost model
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf);
%%% wrap up the constraints as a whole function handle
% path constraint function
path_constraint_func = @(x,u)([u(1)-params.umax;
                               params.umin- u(1);
                               x(2)-1.5*pi;
                               -1.5*pi-x(2);
                               x(1)-0.7;
                               -0.7-x(1)]); % <= 0
Npath_Nineq = numel(path_constraint_func(params.x0, zeros(params.nu,1)));

% bnd constraint function
final_constraint_func = @(xf)([xf(1)-params.xf(1);
                               params.xf(1)-xf(1);
                               xf(2)-params.xf(2);
                               params.xf(2)-xf(2)]); % < 0
Nfinal_Nineq = numel(final_constraint_func(params.xf));

% cost function
obj_path_func = @(x, u, xref, uref) (1/2*(x-xref).'*params.Q*(x-xref) + 1/2*(u-uref).'*params.R*(u-uref)) * params.dt;
obj_bnd_func = @(x) 1/2*(x-params.xf).'*params.Qf*(x-params.xf);                                

% create constraint objects
path_cons = path_constraint(Npath_Nineq);
final_cons = final_constraint(Nfinal_Nineq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[success] = Setup_Functions_barrier(rbt,...
                                    path_constraint_func, final_constraint_func, ...
                                    obj_path_func, obj_bnd_func, ...
                                    params);

%% ................................................................
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

%% Solve ................................
tstart = tic;
[xsol, usol, Ksol] = solver.Solve(rbt,cost,path_cons,final_cons,params);
telapsed = toc(tstart)

J_hist = solver.Jstore;
Jr_hist = solver.J_real;
R_hist = solver.Contract_Rate;
real_cost = compute_cost(xsol,usol,cost,params);
%% plot
figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0,'MarkerSize',3);
ha=gca;
set(ha,'yscale','log');
set(ha,'xscale','log');
grid on;

figure(8888);
plot(Jr_hist,'r-o','LineWidth',2.0,'MarkerSize',3);
ha=gca;
set(ha,'yscale','log');
set(ha,'xscale','log');
grid on;

figure(999);
subplot(2,1,1);
plot(t, (xsol(1:2,:)),'LineWidth',2.0);hold on;
legend('x1', 'x2','Interpreter','latex','FontSize',15);
grid on;

subplot(2,1,2);
plot(t, (xsol(3:4,:)),'LineWidth',2.0);hold on;
legend('dx1', 'dx2','Interpreter','latex','FontSize',15);
grid on;

figure(1000);
yyaxis left
p1=plot(t(1:end-1), usol(1,:),'Color',[0.8 0 0.0],'LineWidth',2.0);hold off;
ylabel('$u_1$','Interpreter','latex','FontSize',15);
yyaxis right
for kk=1:size(Ksol,2)  
    p = plot(t(1:end-1), squeeze(Ksol(1,kk,:)),'Color','r','LineWidth',2.0);hold on;
    p.Color(4)=0.1;
    set(p,'LineStyle','-');
end
ylabel('$u_1$ Gains','Interpreter','latex','FontSize',15);
grid on;

constraint_voilation = zeros(Npath_Nineq, params.N);
for m=1:params.N
    xm = xsol(:,m);
    um = usol(:,m);
    constraint_voilation(:,m) = path_constraint_func(xm,um);
end
figure(3000);
plot(t(1:end-1),constraint_voilation,'LineWidth',2.0); hold on;
plot(t(1:end-1),0*constraint_voilation,'k-.','LineWidth',2.0);
grid on;
%%
%%%% data logging %%%
if logger == 1
    file_name1 = strcat('.\data\T_', ... 
                         num2str(params.shooting_phase),'_', exp_date);
    file_name2 = strcat('.\data\M_', ...
                         num2str(params.shooting_phase), '_', exp_date);
    file_name3 = strcat('.\data\R_', ...
                         num2str(params.shooting_phase), '_', exp_date);
    save(file_name1,'telapsed');
    save(file_name2,'J_hist');
    save(file_name3,'R_hist');
end
