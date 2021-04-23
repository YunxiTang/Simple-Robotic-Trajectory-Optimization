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
params.dt    = .01;
params.T     = 5.0;
params.N     = params.T / params.dt;
params.shooting_phase = 50;
params.x0    = [-0.5;0.0;0.0;0.0];
params.xf    = [2.0;2.0;pi/2;0.0];
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([1 1 1 1])*2;
params.R     = diag([0.1 0.1]);
params.Qf    = diag([5 5 5 5])*8;
params.Rf    = eye(params.nu);
params.Reg_Type = 2;        % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 4.5;
params.umin  = -4.5;
params.Debug = 1;           % 1 -> show details
params.plot = 1;            % 1 -> show plots during optimization
params.Max_iter = 1000;
params.stop = 1e-7;
params.qp = 1;
params.warm_start = 1;
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
car = Car();

%%% wrap up the constraints as a whole function handle
% path constraint
Obstacles = [0.0 1.3 2.0;
             1.0 1.3 2.5;
             0.5 0.6 0.4];

path_constraint_func = @(x,u)([ 0.5*0.5 - ((x(1)-0.0)*(x(1)-0.0) + (x(2)-1.0)*(x(2)-1.0));
                                0.6*0.6 - ((x(1)-1.3)*(x(1)-1.3) + (x(2)-1.3)*(x(2)-1.3));
                                0.4*0.4 - ((x(1)-2.0)*(x(1)-2.0) + (x(2)-2.5)*(x(2)-2.5));
                                u(1)-4.5;
                               -u(1)-4.5;
                                u(2)-4.5;
                               -u(2)-4.5]); % <= 0

% TO DO: add final state constraint here
final_constraint_func = @(xf)([xf(1)-params.xf(1);
                               params.xf(1)-xf(1);
                               xf(2)-params.xf(2);
                               params.xf(2)-xf(2)]); % < 0

% cost model
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions_barrier(params, car, cost, path_constraint_func, final_constraint_func);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

% create constraint objects
path_cons = path_constraint();
final_cons = final_constraint();

%% Solve ................................
tstart = tic;
[xsol, usol, Ksol] = solver.Solve(car,cost,path_cons,final_cons,params);
telapsed = toc(tstart);

J_hist = solver.Jstore;
R_hist = solver.Contract_Rate;
%% plot
figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0,'MarkerSize',3);
ha=gca;
set(ha,'yscale','log');
set(ha,'xscale','log');
grid on;

figure(999);
plot(t, (xsol),'LineWidth',2.0);hold on;
legend("$x$","$y$","$\theta$","$v$",'Interpreter','latex','FontSize',12);
grid on;

figure();
plot(t(1:end-1), usol,'LineWidth',2.0);
legend("$u_{\theta}$","$u_v$",'Interpreter','latex','FontSize',12);
grid on;

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


figure(2000);
plot_obstacle(Obstacles, 2000);hold on;
plot(params.x0(1), params.x0(2), 'kp', 'MarkerFaceColor', 'b', 'MarkerSize', 15); hold on;
plot(params.xf(1), params.xf(2), 'rh', 'MarkerFaceColor', 'r', 'MarkerSize', 15); hold on;
plot(xsol(1,:),xsol(2,:),'r-.','LineWidth',2.0);
h=legend('$Obstacle\;1$','$Obstacle\;2$','$Obstacle\;3$','$Start \; Point$','$Goal\;Point$','$Trajectory$', 'Interpreter','latex','FontSize',13);
h.NumColumns = 2;
set (gcf,'Position',[400,100,500,500], 'color','w');
% grid on;

constraint_voilation = zeros(7, params.N);
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
