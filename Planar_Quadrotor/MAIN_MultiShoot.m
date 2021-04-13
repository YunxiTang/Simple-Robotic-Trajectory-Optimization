%%% Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%    Date  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_date = '23-Mar-2021';
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
log = 0;
params.dt    = .01;
params.T     = 6.0;
params.N     = params.T / params.dt;
params.shooting_phase = 20;
params.x0    = [4.0; 0.0; -0.8; 0.0; 0.0; 0.0];
params.xf    = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0];
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([0.1 0.1 0.1 0.1 0.1 0.1])*5;
params.R     = diag([0.1 0.1])*20;
params.Qf    = diag([50 50 50 50 50 50]);
params.Rf    = eye(params.nu);
params.Reg_Type = 1.0;                    % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 4.4;
params.umin  = 1.0;
params.Debug = ~log;     % 1 -> show details
params.plot = ~log;      % 1 -> show plots during optimization
params.Max_iter = 500;
params.stop = 1e-9;
params.qp = 1;        % 1 -> BoxQP for input constraint
params.clamp = 0;        % 1 -> clamp for input constraint
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
subplot(2,1,1);
yyaxis left
p1=plot(t(1:end-1), usol(1,:),'Color',[0.8 0 0.0],'LineWidth',2.0);hold off;
ylabel('$u_1$','Interpreter','latex','FontSize',15);
yyaxis right
for kk=1:size(Ksol,2)-2  
    p = plot(t(1:end-1), squeeze(Ksol(1,kk,:)),'Color','r','LineWidth',2.0);hold on;
    p.Color(4)=0.1;
    set(p,'LineStyle','-');
end
ylabel('$u_1$ Gains','Interpreter','latex','FontSize',15);
grid on;

subplot(2,1,2);
yyaxis left
plot(t(1:end-1), usol(2,:),'Color',[0 0 0.8],'LineWidth',2.0);hold off;
ylabel('$u_2$','Interpreter','latex','FontSize',15);
for kk=1:size(Ksol,2)-2
    yyaxis right
    pkk = plot(t(1:end-1), squeeze(Ksol(2,kk,:)),'LineStyle','-','Color','b','LineWidth',2.0);hold on;
    pkk.Color(4)=0.1;
end
ylabel('$u_2$ Gains','Interpreter','latex','FontSize',15);
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