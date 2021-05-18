%%% Multiple shooting SLQ for free robot dynamics
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
params.shooting_phase   = 1000;
params.L                = (params.N - 1) / params.shooting_phase + 1;
params.x0               = [0.0;3.5;1.9;0.0;0.2;0.0];
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
params.Q                = diag([120 120 10 50 50 1]);
params.R                = diag([1 1]);
params.Qf               = diag([120 120 10 50 50 1]);
params.Rf               = eye(params.nu);
params.Reg_Type         = 2;     % 1->reg of Quu  / 2->reg of Vxx
params.umax             = 15.0;
params.umin             = 0.01;
params.Debug            = 1;     % 1 -> show details
params.plot             = 0;     % 1 -> show plots during optimization
params.Max_iter         = 1000;
params.stop             = 1e-5;
params.qp               = 1;
params.clamp            = 0;
nt                      = params.T / params.shooting_phase;
tax                     = cell(params.shooting_phase,1);
for i=1:params.shooting_phase
    tax{i}=linspace((i-1)*nt,i*nt,nt/params.dt+1);
end
params.t = tax;

params.MapNo = 1;
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
telapsed = toc(tstart);
%%
figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
J_hist = solver.Jstore;

figure(666);
plot(tref, xsol(3,:),'LineWidth',2.0);

figure(999);
plot(xsol(1,:), xsol(2,:),'r--','LineWidth',2.0); hold on;
plot(xref(1,:), xref(2,:), 'k-', 'LineWidth',2.0);
grid on;

figure(1000);
subplot(2,1,1);
yyaxis left
p1=plot(tref(1:end-1), usol(1,:),'Color',[0.8 0 0.0],'LineWidth',2.0);hold off;
ylabel('$u_1$','Interpreter','latex','FontSize',15);
yyaxis right
% for kk=1:size(Ksol,2)-2  
%     p = plot(t(1:end-1), squeeze(Ksol(1,kk,:)),'Color','r','LineWidth',2.0);hold on;
%     p.Color(4)=0.1;
%     set(p,'LineStyle','-');
% end
ylabel('$u_1$ Gains','Interpreter','latex','FontSize',15);
grid on;

subplot(2,1,2);
yyaxis left
plot(tref(1:end-1), usol(2,:),'Color',[0 0 0.8],'LineWidth',2.0);hold off;
ylabel('$u_2$','Interpreter','latex','FontSize',15);
% for kk=1:size(Ksol,2)-2
%     yyaxis right
%     pkk = plot(t(1:end-1), squeeze(Ksol(2,kk,:)),'LineStyle','-','Color','b','LineWidth',2.0);hold on;
%     pkk.Color(4)=0.1;
% end
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
%%
%%%%%%%%% animation %%%%%%%%%
figure(2000);
k = 25;
plot(xsol(1,:),xsol(2,:),'r-.','LineWidth',2.0);
hold off;
planar_quad.animation(tref,xsol,k,2000);
set (gcf,'Position',[400,100,500,500], 'color','w');
axis equal;