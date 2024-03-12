function [Xtraj,Utraj,x_dft,J] = forward_multi_rollout1(x_init,x_goal,model,cost,...
                                                        u_max,M,S,dt,TRAIL)
%FORWARD_MULTI_ROLLOUT forward rollout for multiple shooting.
%%%        ARGS                                     DIMENSIONS
%%% INPUT: x_init:          initial x               [Nx, 1]
%%%        x_goal:          goal x                  [Nx, 1] 
%%%        model:           robot model             [class]
%%%        cost:            cost model              [class]
%%%        umax:            max torque              [Nu, 1]
%%%        M:               number of control inputs    [1]
%%%        S:               number of segments          [1]
%%%        dt:              grid time interval          [1]
%%%        TRAIL:           sol_mdl                 [class]
Nx = numel(x_init);
Nu = numel(u_max);
L = M / S;
Xtraj = zeros(S, L+1 , Nx, 1);
Utraj = 5.0 * ones(S, L   , Nu, 1);
T = zeros(S, L+1, 1);
% first forward simulation
x_S_guess = [linspace(x_init(1),x_goal(1),S+1);
             linspace(x_init(2),x_goal(2),S+1);
             linspace(x_init(3),x_goal(3),S+1);
             linspace(x_init(4),x_goal(4),S+1)];

% for init value for each multiple shooting interval.
t_now = 0.0;
k = 1;
for i=1:S
    Xtraj(i,1,:,:) = x_S_guess(:,i);
    T(i,1,:) = t_now;
    for j=1:L
        t_now = k * dt;
        T(i,j+1,:) = t_now;
        % rollout with rk45
        Xtraj(i,j+1,:,:) = model.rk45(Xtraj(i,j,:,:), Utraj(i,j,:,:), dt);
        k = k + 1;
    end
end
TRAIL.store_time(T);
TRAIL.store_trial(Xtraj,Utraj);
x_dft = TRAIL.store_defect(Xtraj);
plot_trail(TRAIL.time_bar,TRAIL.XTrials{1},TRAIL.Defects{1},'$Init\;Trail$');
J = cost.calc_J(Xtraj,Utraj,x_goal,zeros(Nu,1),dt);
end

