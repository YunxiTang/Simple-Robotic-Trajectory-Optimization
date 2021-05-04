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
params.dt    = .01;
params.T     = 2.0;
params.N     = params.T / params.dt;
params.shooting_phase = 1;
params.x0    = [0.0; 0.0; 0.0; 0.0];
params.xf    = [pi;  0.0; 0.0; 0.0];
params.nx    = numel(params.x0);
params.nu    = 2 ;
params.Q     = diag([0.1 0.1 0.1 0.1])*10;
params.R     = 1* eye(params.nu);
params.Qf    = diag([0.1 0.1 0.1 0.1])*1200;
params.Rf    = eye(params.nu);
params.Reg_Type = 2;  % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 7;
params.umin  = -7;
params.Debug = 1;     % 1 -> show details
params.plot = 1;      % 1 -> show plots during optimization
params.Max_iter = 2000;
params.stop = 1e-7;
params.qp = 1;        % 1 -> BoxQP for input constraint
params.clamp = 0;        % 1 -> clamp for input constraint
nt = params.T / params.shooting_phase;
tax = cell(params.shooting_phase,1);
for i=1:params.shooting_phase
    tax{i}=linspace((i-1)*nt,i*nt,nt/params.dt+1);
end
params.t = tax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rr_bt = RR_robot();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, rr_bt, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

tstart = tic;
[xsol, usol, Ksol] = solver.Solve(rr_bt,cost,params);
telapsed = toc(tstart)

figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
J_hist = solver.Jstore;
R_hist = solver.Contract_Rate;
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
t = 0.0:params.dt:params.T;
figure(789);
plot(t,xsol,'LineWidth',2.0);
legend("$x_1$","$x_2$","$x_3$","$x_4$",'Interpreter','latex','FontSize',12);
figure(799);
plot(t(1:end-1),usol,'LineWidth',2.0);
legend("$u_1$","$u_2$",'Interpreter','latex','FontSize',12);

figure(1000);
subplot(2,1,1);
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

subplot(2,1,2);
yyaxis left
plot(t(1:end-1), usol(2,:),'Color',[0 0 0.8],'LineWidth',2.0);hold off;
ylabel('$u_2$','Interpreter','latex','FontSize',15);
for kk=1:size(Ksol,2)
    yyaxis right
    pkk = plot(t(1:end-1), squeeze(Ksol(2,kk,:)),'LineStyle','-','Color','b','LineWidth',2.0);hold on;
    pkk.Color(4)=0.1;
end
ylabel('$u_2$ Gains','Interpreter','latex','FontSize',15);
grid on;
% figure(790);
% plot(t(1:end-1),squeeze(Ksol(1,:,:)),'LineWidth',2.0);
% legend("$K_1$","$K_2$","$K_3$","$K_4$",'Interpreter','latex','FontSize',12);