%%% Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    Date  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_date = '01-Apr-2021';
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
log = 0;
params.dt    = .01;
params.T     = 10.0;
params.N     = params.T / params.dt;
params.shooting_phase = 10;
params.x0    = [0.5;0.5;0.5;0.0;0.0;0.00;0.0;0.0;0.0;0.0];
params.xf    = zeros(10,1);
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([1;1;1;0.00;0.00;1;1;1;0.00;0.00])*20;
params.R     = eye(params.nu)*0.1;
params.Qf    = diag([1;1;1;0.0;0.0;1;1;1;0.0;0.0])*500;
params.Rf    = eye(params.nu);
params.Reg_Type = 2;     % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 0.4;
params.umin  = -0.4;
params.Debug = 1;     % 1 -> show details
params.plot = 1;      % 1 -> show plots during optimization
params.Max_iter = 100;
params.stop = 1e-7;
params.qp = 1;
nt = params.T / params.shooting_phase;
tax = cell(params.shooting_phase,1);
for i=1:params.shooting_phase
    tax{i}=linspace((i-1)*nt,i*nt,nt/params.dt+1);
end
params.t = tax;
t = 0.0:params.dt:params.T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fallcat = falling_cat();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup_Functions(params, fallcat, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

tstart = tic;
[xsol, usol, Ksol] = solver.Solve(fallcat,cost,params);
telapsed = toc(tstart)

figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
J_hist = solver.Jstore;
R_hist = solver.Contract_Rate;

figure(999);
subplot(2,1,1);
plot(t, rad2deg(xsol(1:3,:)),'LineWidth',2.0);hold on;
plot(t, rad2deg(xsol(6:8,:)),'LineStyle','-.','LineWidth',2.0);
legend('Body Yaw', 'Body Pitch', 'Body Roll','Interpreter','latex','FontSize',15);
grid on;
subplot(2,1,2);
plot(t, rad2deg(xsol(4:5,:)),'LineWidth',2.0);hold on;
plot(t, rad2deg(xsol(9:10,:)),'LineStyle','-.','LineWidth',2.0);
legend('Tail Yaw', 'Tail Pitch','Interpreter','latex','FontSize',15);
grid on;

figure(1000);
yyaxis left;
plot(t(1:end-1), usol(1,:),'r','LineWidth',2.0);hold on;
plot(t(1:end-1), usol(2,:),'b','LineWidth',2.0);
% legend('$u_1$','$u_2$','Interpreter','latex','FontSize',15);
yyaxis right;
plot(t(1:end-1), squeeze(Ksol(1,:,:)),'r-.','LineWidth',2.0);hold on;
plot(t(1:end-1), squeeze(Ksol(2,:,:)),'b-.','LineWidth',2.0);hold on;
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
