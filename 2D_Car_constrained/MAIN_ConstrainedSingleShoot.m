%%% Constrained Single shooting SLQ for free robot dynamics
%%% Y.X TANG (yxtang@mae.cuhk.edu.hk BMT LAB, CUHK)
clc;
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
log = 1;       % data log flag  
params.dt    = .05;
params.T     = 5.0;
params.N     = params.T / params.dt;
params.x0    = [0.0;0.0;0.0;0.0];
params.xf    = [3.0;3.0;-pi/2;0.0];
params.nx    = numel(params.x0);
params.nu    = 2;
params.Q     = diag([0.1 0.1 0.1 1.0]);
params.R     =  diag([0.2 0.1]);
params.Qf    =  diag([50 50 50 50]);
params.Rf    = eye(params.nu);
params.Reg_Type = 2;  % 1->reg of Quu  / 2->reg of Vxx
params.umax  = 5.0;
params.umin  = -5.0;
params.Debug = 1;     % 1 -> show details
params.plot = 1;      % 1 -> show plots during optimization
params.Max_iter = 500;
params.stop = 1e-9;
taxis = linspace(0,params.T,params.N);
params.tax = taxis;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create Obstacles on the Map %%%
%%%% --------------------------  %%%
%%%% |xc_1 | xc_2 | ... | xc_m|  %%%
%%%% |-----+------+-----+-----|  %%%
%%%% |yc_1 | yc_2 | ... | yc_m|  %%%
%%%% |-----+------+-----+-----|  %%%
%%%% | r1  |  r2  | ... |  rm |  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Obstacles = [2.0  0.5  3.0;
             2.0  1.0  2.4;
             0.4  0.3  0.15];
Constraints = constraint(Obstacles);
% show map
Map = 1;
figure(Map);
plot(params.x0(1), params.x0(2), 'kp', 'MarkerFaceColor', 'b', 'MarkerSize', 15); hold on;
plot(params.xf(1), params.xf(2), 'rh', 'MarkerFaceColor', 'r', 'MarkerSize', 15); hold on;
Constraints.plot_obstacle(Map);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% create robot and cost mdl %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
car = Car();
cost = cst_mdl(params.Q,params.R,params.Qf,params.Rf,params.umax,params.umin);