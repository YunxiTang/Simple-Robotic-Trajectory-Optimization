%%% Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
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
log = 0;
params.dt    = .001;
params.T     = 0.3;
params.N     = params.T / params.dt;
params.shooting_phase = 1;
params.x0    = [0.3;0.3;0.2;0.0;0.0;0.00;0.0;0.0;0.0;0.0];
params.xf    = zeros(10,1);
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([1;1;1;0.00;0.00;1;1;1;0.00;0.00])*5;
params.R     = eye(params.nu)*0.001;
params.Qf    = diag([1;1;1;0.00;0.00;1;1;1;0.00;0.00])*6;
params.Rf    = eye(params.nu);
params.Reg_Type = 1;        % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 0.5;
params.umin  = -0.5;
params.Debug = 1;           % 1 -> show details
params.plot = 1;            % 1 -> show plots during optimization
params.Max_iter = 500;
params.stop = 1e-7;
params.qp = 0;
nt = params.T / params.shooting_phase;
tax = cell(params.shooting_phase,1);
for i=1:params.shooting_phase
    tax{i}=linspace((i-1)*nt,i*nt,nt/params.dt+1);
end
params.t = tax;
t = 0.0:params.dt:params.T;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot,constraints & cost mdl %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% build a robot
fallcat = falling_cat();

%%% wrap up the constraints as a whole function handle
% path constraint
path_constraint_func = @(x,u)([ u(1)-0.5;
                               -u(1)-0.5;
                                u(2)-0.5;
                               -u(2)-0.5;
                                x(4)-deg2rad(50);
                               -x(4)-deg2rad(50);
                                x(5)-deg2rad(40);
                               -x(5)-deg2rad(40)]); % <= 0

% TO DO: add final state constraint here
final_constraint_func = @(xf)([xf(1)-deg2rad(0);
                              -xf(1)-deg2rad(0);
                               xf(2)-deg2rad(0);
                              -xf(2)-deg2rad(0);
                               xf(3)-deg2rad(0);
                              -xf(3)-deg2rad(0)]); % < 0

% cost model
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions_barrier(params, fallcat, cost, path_constraint_func, final_constraint_func);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

% create constraint objects
path_cons = path_constraint();
final_cons = final_constraint();

%% Solve ................................
tstart = tic;
[xsol, usol, Ksol] = solver.Solve(fallcat,cost,path_cons,final_cons,params);
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
legend('Body Yaw', 'Body Pitch', 'Body Roll','Interpreter','latex','FontSize',15);
grid on;

subplot(2,1,2);
plot(t, rad2deg(xsol(4:5,:)),'LineWidth',2.0);hold on;
% plot(t, rad2deg(xsol(9:10,:)),'LineStyle','-.','LineWidth',2.0);
plot(t, 50*ones(params.N+1),'k--','LineWidth',2.0);hold on;
plot(t, -50*ones(params.N+1),'k--','LineWidth',2.0);hold on;
plot(t, 40*ones(params.N+1),'r--','LineWidth',2.0);hold on;
plot(t, -40*ones(params.N+1),'r--','LineWidth',2.0);hold on;
legend('Tail Yaw', 'Tail Pitch','Interpreter','latex','FontSize',15);
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
