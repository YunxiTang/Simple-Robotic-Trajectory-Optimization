function [Xtraj,Utraj,x_dft,J] = forward_multi_rollout2(x_init,x_goal,x_traj,u_traj,x_dft,model,cost,...
                                                        u_max,M,S,dt,alpha,Fx,Fu,Kfb,kff,TRAIL)
%FORWARD_MULTI_ROLLOUT2 forward rollout for multiple shooting.
%%%        ARGS                                     DIMENSIONS
%%% INPUT: x_init:          initial x               [Nx, 1]
%%%        x_goal:          goal x                  [Nx, 1]
%%%        x_traj           state traj(last iter)   [S, L+1, Nx, 1]
%%%        u_traj           control traj(last iter) [S, L  , Nu, 1]
%%%        x_dft            defects along traj      [N, Nx, 1]
%%%        model:           robot model             [class]
%%%        cost:            cost model              [class]
%%%        umax:            max torque              [Nu, 1]
%%%        M:               number of control inputs    [1]
%%%        S:               number of segments          [1]
%%%        dt:              grid time interval          [1]
%%%        alpha:           line search factor          [1]
%%%        TRAIL:           sol_mdl                 [class]
Nx = numel(x_init);
Nu = numel(u_max);
L = M / S;
Xtraj = zeros(S, L+1 , Nx, 1);
Utraj = zeros(S, L   , Nu, 1);
Xtraj(1,1,:,:) = x_init;
% linear sweep to update intermediate node state decision variable
for i=1:S
    % defect at node
    if 1 < i 
        defect = x_dft(i-1,:,:)';
        fx = squeeze(Fx(i-1,end,:,:));
        fx = reshape(fx, [Nx, Nx]);
        fu = squeeze(Fu(i-1,end,:,:));
        fu = reshape(fu, [Nx, Nu]);
        
        x_last_pre_seg = squeeze(x_traj(i-1,end,:,:));
        x_last_pre_seg = reshape(x_last_pre_seg, [Nx, 1]);
        
        x_now_pre_seg = squeeze(Xtraj(i-1,end,:,:));
        x_now_pre_seg = reshape(x_now_pre_seg, [Nx, 1]);
        dx_hat = x_now_pre_seg - x_last_pre_seg; 
%         Xtraj(i,1,:,:) = x_now_pre_seg + 1.0 * (fx * dx_hat + fu * du_star) + 2.0 * defect;
        Xtraj(i,1,:,:) = reshape(squeeze(Xtraj(i-1,end,:,:)),[Nx,1]) - (1 - 0.5) * defect;
    end
    for j=1:L
        x_now = squeeze(Xtraj(i,j,:,:));
        x_now = reshape(x_now, [Nx,1]);
        
        x_last = squeeze(x_traj(i,j,:,:));
        x_last  = reshape(x_last,[Nx, 1]);
        
        u_last = squeeze(u_traj(i,j,:,:));
        u_last = reshape(u_last, [Nu,1]);
        l_n = squeeze(kff(i,j,:,:));
        L_n = squeeze(Kfb(i,j,:,:));
        l_n = reshape(l_n,[Nu, 1]);
        L_n = reshape(L_n,[Nu, Nx]);
        du_ff = alpha / 1.1 * l_n;
        du_fb = L_n * (x_now-x_last);
        
        u_ctr = u_last + du_ff + du_fb;
        Utraj(i,j,:,:) = u_ctr;
        % forward rollout
        Xtraj(i,j+1,:,:) = model.rk45(x_now, u_ctr, dt);
    end
    du_star = du_ff + du_fb;
end
TRAIL.store_trial(Xtraj,Utraj);
J = cost.calc_J(Xtraj,Utraj,x_goal,zeros(Nu,1),dt);
x_dft = TRAIL.store_defect(Xtraj);
end

