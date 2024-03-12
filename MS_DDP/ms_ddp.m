% multiple shooting for ddp solver for optimal contol
function [SOL,K_f,x_rep,JJ] = ms_ddp(x_init,x_goal,model,cost,u_max,N,M,S,dt,...
                                  num_iter,alpha,stop_criterion,TRAIL)
%%%        ARGS                                    DIMENSIONS
%%% INPUT: x_init:               initial state              [nx, 1]
%%%        x_goal:               goal state                 [nx, 1]
%%%        model:                robot model                [class]
%%%        cost:                 cost model                 [class]
%%%        u_max:                torque limitation          [nu,1]
%%%        N:                    number of states           [1]
%%%        M:                    number of control inputs   [1]
%%%        S:                    number of segments         [1]
%%%        dt:                   grid time interval         [1]
%%%        num_iter:             maximum iteration times    [1]
%%%        alpha:                factor for line search     [1]
%%%        stop_criterion:       stop criterion             [1]
%%%        TRAIL:                sol_mdl                    [class]
    %% allocate arrays memory for DDP structure
    persistent L Nx Nu
    if isempty(L) || isempty(Nx) || isempty(Nu)
        L = M / S;
        Nx = numel(x_init);
        Nu = numel(u_max);
    end

    % cost history
    JJ = zeros(num_iter,1);
    
    % state history
    x_traj = zeros(S, L+1, Nx, 1);
    u_traj = zeros(S, L, Nu, 1);
    %% initialize MS_SQL solver with initial control guessed input and state
    % First multiple shooting forward simulation
    [x_traj_new,u_traj_new,x_dft,J] = forward_multi_rollout1(x_init,x_goal,model,cost,...
                                                             u_max,M,S,dt,TRAIL);

    %% perform main MS_SLQ iterations on the trajectory and input sequence
    disp("[INFO]: Executing SLQ Iteration.");
    for i=1:num_iter
        cost.update_iter();
        
        % FORWARD SIMULATION WITH updated control sequence from previous iteration
        if i > 1
            
            [x_traj_new,u_traj_new,x_dft,J] = forward_multi_rollout2(x_init,x_goal,x_traj,u_traj,x_dft,model,cost,...
                                                                     u_max,M,S,dt,alpha,F_x,F_u,Kfb,kff,TRAIL);
            if mod(i, 20) == 0
                plot_trail(TRAIL.time_bar,TRAIL.XTrials{i},TRAIL.Defects{i},'$Trail$');
%                 pause;
            end
        end
        
        % update trajectory
        x_traj = x_traj_new;
        u_traj = u_traj_new;
        
        % compute total cost
        JJ(i) = J;
        fprintf(' -------DDP Iteration %d out of %d  ||  Cost: %.8f------- \n',i,num_iter, JJ(i));
        % convert state trajectory and control trajectory into a whole
        % trail
        flat_x_traj = TRAIL.flatted_traj(x_traj);
        flat_u_traj = TRAIL.flatted_u(u_traj);
        
        % Next Step: perform backwards propogation
        [V,V_x,V_xx,Q_x,Q_u,Q_uu,Q_xx,Q_xu,Q_ux,F_x,F_u,kff,Kfb] = ...
        backward_multi_prop(model, cost, x_traj, u_traj, x_goal, 0, x_dft, dt, 1);
        % when to stop
        if i > 1
            if abs(JJ(i-1)-JJ(i)) <= stop_criterion 
                fprintf('[INFO]: Cost Converge at %d-th iteration with Cost of %7f.\n',i,JJ(i));
                break
            end
        end
    end
    SOL = u_traj;
    K_f = Kfb;
    x_rep = x_traj;
end
