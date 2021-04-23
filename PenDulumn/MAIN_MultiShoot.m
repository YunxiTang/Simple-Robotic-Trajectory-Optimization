%%% Multiple shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
log = 1;
params.dt    = .01;
params.T     = 5.0;
params.N     = params.T / params.dt;
params.shooting_phase = 5;
params.x0    = [0.0;0.0];
params.xf    = [3.14;0.0];
params.nx    = numel(params.x0);
params.nu    = 1;
params.Q     = eye(params.nx);
params.R     = 1*eye(params.nu);
params.Qf    = 100*eye(params.nx);
params.Rf    = eye(params.nu);
params.Reg_Type = 1;  % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 3.5;
params.umin  = -3.5;
params.Debug = 0;     % 1 -> show details
params.plot = 0;      % 1 -> show plots during optimization
params.Max_iter = 1000;
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
[xbar,ubar,du,K,dft, xsol, usol, Ksol] = solver.Solve(pendulum,cost,params);
telapsed = toc(tstart)

figure(888);
plot(solver.Jstore,'b-o','LineWidth',2.0);
J_hist = solver.Jstore;
if log == 1
    file_name1 = strcat('.\data\pd_msddp\T',num2str(params.shooting_phase));
    file_name2 = strcat('.\data\pd_msddp\M',num2str(params.shooting_phase));
    save(file_name1,'telapsed');
    save(file_name2,'J_hist');
end

if log == 2
    file_name3 = strcat('.\model_sensitve_analysis\usol',num2str(params.shooting_phase));
    file_name4 = strcat('.\model_sensitve_analysis\xsol',num2str(params.shooting_phase));
    save(file_name3,'usol');
    save(file_name4,'xsol');
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