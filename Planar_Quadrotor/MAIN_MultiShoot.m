%%% Multiple shooting SLQ for free robot dynamics
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
params.dt    = .02;
params.T     = 4.0;
params.N     = params.T / params.dt;
params.shooting_phase = 20;
params.x0    = [1.0; 0.2; 1.0; 0.0; 0.0; 0.0];
params.xf    = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([0.1 0.1 0.1 0.1 0.1 0.1]);
params.R     = diag([0.1 0.1]);
params.Qf    = diag([50 50 50 50 50 50])/5;
params.Rf    = eye(params.nu);
params.Reg_Type = 1.0;                    % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 5;
params.umin  = -5;
params.Debug = 1;     % 1 -> show details
params.plot = 1;      % 1 -> show plots during optimization
params.Max_iter = 500;
params.stop = 1e-9;
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
planar_quad = planar_quadrotor();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, planar_quad, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

tstart = tic;
[xsol, usol, Ksol] = solver.Solve(planar_quad,cost,params);
telapsed = toc(tstart)

figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
J_hist = solver.Jstore;

figure(999);
plot(t, xsol,'LineWidth',2.0);
grid on;

figure(1000);
plot(t(1:end-1), usol,'LineWidth',2.0);
grid on;

%%%% data logging %%%
if log == 1
    file_name1 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\Planar_Quadrotor\data\T_', ... 
                         num2str(params.shooting_phase),'_', exp_date);
    file_name2 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\Planar_Quadrotor\data\M_', ...
                         num2str(params.shooting_phase), '_', exp_date);
    save(file_name1,'telapsed');
    save(file_name2,'J_hist');
end

%%%%%%%%% animation %%%%%%%%%
% figure(789);
% k = 2;
% plot(xsol(1,:),xsol(2,:),'k-.','LineWidth',1.0); hold on;
% planar_quad.animation(t,xsol,k,789);
% plot(xsol(1,:),xsol(2,:),'k-.','LineWidth',1.0); hold on;
% plot(xsol(1,1:k:end),xsol(2,1:k:end),'ro','LineWidth',2.0); hold on;