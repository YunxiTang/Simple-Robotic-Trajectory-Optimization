%%% Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
log = 0;
params.dt    = .01;
params.T     = 5.0;
params.N     = params.T / params.dt;
params.shooting_phase = 10;
params.x0    = [0.0;0.0];
params.xf    = [3.14;0.0];
params.nx    = numel(params.x0);
params.nu    = 1;
params.Q     = eye(params.nx);
params.R     = eye(params.nu);
params.Qf    = 10*eye(params.nx);
params.Rf    = eye(params.nu);
params.Reg_Type = 2;  % 1->reg of Quu  / 2->reg of Vxx
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pendulum = rbt_mdl();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, pendulum, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

tstart = tic;
[xbar,ubar,du,K,dft] = solver.Solve(pendulum,cost,params);
telapsed = toc(tstart)

figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
J_hist = solver.Jstore;
if log == 1
    file_name1 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\MIXSOL\data\msddp\T',num2str(params.shooting_phase));
    file_name2 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\MIXSOL\data\msddp\M',num2str(params.shooting_phase));
    save(file_name1,'telapsed');
    save(file_name2,'J_hist');
end

for k=1:params.shooting_phase
    if mod(k,2)==0
        clr = 'r';
    else
        clr = 'b';
    end
    figure(345);hold on;
    plot(tax{k},ubar{k},'Color',clr,'LineWidth',2.0);
    xlabel('Time[s]', 'Interpreter','latex','FontSize',15);
    ylabel('Torque[N]','Interpreter','latex','FontSize',15);
    grid on;
    title('Control','Interpreter','latex','FontSize',15);
    figure(456);hold on;
    plot(tax{k},xbar{k},'Color',clr,'LineWidth',2.0);
    grid on;
end