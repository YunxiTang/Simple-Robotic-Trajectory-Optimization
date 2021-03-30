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
log = 1;
params.dt    = .02;
params.T     = 5.0;
params.N     = params.T / params.dt;
params.shooting_phase = 5;
params.x0    = [-2.0;2.0;0.0;0.0];
params.xf    = [4.0;4.0;0.0;0.0];
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([.5 0.5 0.1 0.1]);
% params.R     = 0.1 * eye(params.nu);
% params.Qf    =  50 * eye(params.nx);
params.R     =  diag([0.2 0.1]);
params.Qf    =  diag([20 20 20 20]);
params.Rf    = eye(params.nu);
params.Reg_Type = 2;  % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 5;
params.umin  = -5;
params.Debug = 0;     % 1 -> show details
params.plot = 1;      % 1 -> show plots during optimization
params.Max_iter = 3000;
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
car = Car();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function Setup %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Setup_Functions(params, car, cost);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Call Solver %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver = msddp_solver(params);

tstart = tic;
[xbar,ubar,du,K,dft] = solver.Solve(car,cost,params);
telapsed = toc(tstart)

figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
J_hist = solver.Jstore;
%%%% data logging %%%
if log == 1
    file_name1 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\2D_Car\data\T_', ... 
                         num2str(params.shooting_phase),'_', exp_date);
    file_name2 = strcat('D:\TANG Yunxi\Motion Planning Locomotion\motion_planning\2D_Car\data\M_', ...
                         num2str(params.shooting_phase), '_', exp_date);
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