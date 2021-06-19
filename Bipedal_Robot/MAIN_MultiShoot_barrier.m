%%% Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;
close all;
addpath('.\visualize', '.\kinematics');
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    Date  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% 
log = 0;
params.dt    = .001;
params.T     =  0.2;
params.N     = params.T / params.dt;
params.shooting_phase = 1;
params.x0    = [-deg2rad(4);  deg2rad(4); deg2rad(3); [0.8;-0.8;0.0]];
params.xf    = [ deg2rad(4); -deg2rad(4); deg2rad(0); [0.0; 0.0;0.0]];
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([1;1;0;0;0;0])*1;
params.R     = eye(params.nu)*0.1;
params.Qf    = diag([1;1;0;0;0;0])*10;
params.Rf    = eye(params.nu);
params.Reg_Type = 2;        % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 40;
params.umin  = -40;
params.Debug = 1;           % 1 -> show details
params.plot = 0;            % 1 -> show plots during optimization
params.Max_iter = 1000;
params.stop = 1e-6;
params.qp = 0;
nt = params.T / params.shooting_phase;
tax = cell(params.shooting_phase,1);
for i=1:params.shooting_phase
    tax{i}=linspace((i-1)*nt,i*nt,nt/params.dt+1);
end
params.t = tax;
t = 0.0:params.dt:params.T;
eps = .001;
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot,constraints & cost mdl %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% build a robot
rbt = Bipedal();

%%% wrap up the constraints as a whole function handle
% path constraint
path_constraint_func = @(x,u)([ u(1)-params.umax;
                                params.umin-u(1);
                                u(2)-params.umax;
                                params.umin-u(2);
                                x(2)-(-x(1));
                                (-x(1))-x(2);
                                x(3)-deg2rad(10);
                                deg2rad(0)-x(3);
                                -ground_clear(x, rbt)]); % <= 0

% TO DO: add final state constraint here
final_constraint_func = @(xf)([xf(1)-params.xf(1);
                               params.xf(1)-xf(1);
                               xf(2)-params.xf(2);
                               params.xf(2)-xf(2);
%                                xf(3)-params.xf(3);
%                                params.xf(3)-xf(3);
                               ]); % < 0

% cost model
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions_barrier(params, rbt, cost, path_constraint_func, final_constraint_func);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

% create constraint objects
path_cons = path_constraint();
final_cons = final_constraint();

%% Solve ................................
tstart = tic;
[xsol, usol, Ksol] = solver.Solve(rbt,cost,path_cons,final_cons,params);
telapsed = toc(tstart);
%% plot
figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
ha=gca;
set(ha,'yscale','log');
set(ha,'xscale','log');
grid on;
J_hist = solver.Jstore;
R_hist = solver.Contract_Rate;

figure(999);
subplot(2,1,1);
plot(t, rad2deg(xsol(1:3,:)),'LineWidth',2.0);hold on;
% plot(t, rad2deg(xsol(6:8,:)),'LineStyle','-.','LineWidth',2.0);hold on;
legend('q1', 'q2', 'q3','Interpreter','latex','FontSize',15);
grid on;

subplot(2,1,2);
plot(t, rad2deg(xsol(4:6,:)),'LineWidth',2.0);hold on;
% plot(t, rad2deg(xsol(9:10,:)),'LineStyle','-.','LineWidth',2.0);
legend('dq1', 'dq2', 'dq3', 'Interpreter','latex','FontSize',15);
grid on;

figure(1000);
% yyaxis left;
plot(t(1:end-1), usol(1,:),'r','LineWidth',2.0);hold on;
plot(t(1:end-1), usol(2,:),'b','LineWidth',2.0);hold on;
plot(t(1:end-1), params.umax*ones(params.N,1),'k-.','LineWidth',2.0); hold on;
plot(t(1:end-1), params.umin*ones(params.N,1),'m-.','LineWidth',2.0); hold on;
legend('$u_L$', '$u_R$','$u_{ub}$','$u_{lb}$','Interpreter','latex','FontSize',15);
grid on;
% legend('$u_1$','$u_2$','Interpreter','latex','FontSize',15);
% yyaxis right;
% plot(t(1:end-1), squeeze(Ksol(1,:,:)),'r-.','LineWidth',2.0);hold on;
% plot(t(1:end-1), squeeze(Ksol(2,:,:)),'b-.','LineWidth',2.0);hold on;
% grid on;
%%
N_ineq = numel(path_constraint_func(zeros(params.nx), zeros(params.nu)));
constraint_voilation = zeros(N_ineq, params.N);
for m=1:params.N
    xm = xsol(:,m);
    um = usol(:,m);
    constraint_voilation(:,m) = path_constraint_func(xm,um);
end
figure(3000);
plot(t(1:end-1),constraint_voilation,'LineWidth',2.0); hold on;
plot(t(1:end-1),0*constraint_voilation,'k-.','LineWidth',2.0);
grid on;
%%%% data logging %%%
if log == 1
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

%% 
sln.T = {t};
sln.Y = {xsol};
rbt_animate(sln, rbt)
set (gcf,'Position',[400,100,500,500], 'color','w');
axis off;
title('$Bipedal \; Walking$','Interpreter','latex','FontSize',25);
xlabel('$x$','Interpreter','latex','FontSize',20);
ylabel('$z$','Interpreter','latex','FontSize',20);
%% \
X_swf = zeros(params.N+1,1);
Z_swf = zeros(params.N+1,1);
for i=1:(params.N+1)
    qi = xsol(:,i);
    [x_swf, z_swf, ~, ~] = kin_swf(qi, rbt);
    X_swf(i) = x_swf;
    Z_swf(i) = z_swf;
end

figure();
plot(X_swf,Z_swf);