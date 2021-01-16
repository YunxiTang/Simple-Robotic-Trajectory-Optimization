% ddp solver for fast nonlinear MPC
function [SOL,K_f,x_rep] = ddp_tracking(x_now,x_ref,u_ref,n_h,model,cost,K_lqr,...
                                        u_max,dt,num_iter,alpha,stop_criterion)
%%%        ARGS                                    DIMENSIONS
%%% INPUT: x_now:       current state              [nx, 1]
%%%        x_ref:       ref state trajectory       [nx, n_h]
%%%        u_ref:       ref input trajectory       [nu, n_h]
%%%        n_h:         length of horizon          [1]
%%%        model:       robot model                [class]
%%%        cost:        cost model                 [class]
%%%        K_lqr:       ref feedback policy        [nx, n_h]
%%%        u_max:       torque limitation          [nu,1]
%%%        num_iter:    maximum iteration times    [1]
%%%        alpha:       factor for line search     [1]
%%%        stop_criterion: stop criterion          [1]
    %% allocate arrays memory for DDP structure
    % cost history
    J = zeros(num_iter,1);
    % state history
    x = cell(n_h,1);
    x_new = cell(n_h,1);
    
    x{1} = x_now; 
    x_new{1} = x_now;
    
    % control input his
    u = cell(n_h-1,1);
    
    % Value function and derivatives
    V = zeros(n_h,1); V_x = cell(n_h,1); V_xx = cell(n_h,1);
    
    % State action value derivatives
    Q_x = cell(n_h,1);  Q_u = cell(n_h,1);
    Q_xx = cell(n_h,1); Q_uu = cell(n_h,1);
    Q_xu = cell(n_h,1); Q_ux = cell(n_h,1);
    
    %% initialize ddp solver with local lqr controller
    u{n_h-1} = zeros(numel(u_max),1);
    K{n_h-1} = zeros(numel(x_now),1);
    for i=1:n_h-1
        u{i} = u_ref(i,:) - K_lqr(i,:)*((x_new{i} - x_ref(i,:).'));
        x_new{i+1} = model.rk45(x_new{i},u{i},dt);
    end
    %% perform main DDP iterations on the trajectory and input sequence
    for i=1:num_iter
        cost.update_iter();
        % FORWARD SIMULATION WITH updated control sequence from previous iteration
        if i > 1
            for k=1:n_h-1
                % compute control update feed-forward and feed-back
                % w/ regularzation in case u change a lot
                du_ff = -inv(Q_uu{k}+ 2.5 * eye(numel(u{k}))) * Q_u{k};
                du_fb = -inv(Q_uu{k}+ 2.5 * eye(numel(u{k}))) * Q_ux{k} * (x_new{k} - x{k});
                K{k} = -inv(Q_uu{k}+  2.5 * eye(numel(u{k}))) * Q_ux{k};
                % limit feed-forward control modification with clamping
%                 for m = 1:numel(u_max)
%                     du_ff(m) = min(u_max(m), max(-u_max(m), ...
%                                              du_ff(m) + u{k}(m))) - u{k}(m);
%                 end
                % update control
                u{k} = u{k} + alpha / (1.5) .* (du_ff) + du_fb;
                % FORWARD SIMULATION: Compute next state in trajectory with new control
                x_new{k+1} = model.rk45(x_new{k}, u{k}, dt);
            end 
        end
        
        % update trajectory
        x = x_new;
        
        % compute total cost
        J(i) = cost.phi(x{n_h},x_ref(n_h,:)');
        
        for k=1:n_h-1
            J(i) = J(i) + cost.L(x{k}, u{k}, x_ref(n_h,:)', u_ref(n_h,:)', dt);
        end
        
        % compute terminal value function and derivatives
        V(n_h) = cost.phi(x{n_h},x_ref(n_h,:)');
        V_x{n_h} = cost.phi_x(x{n_h},x_ref(n_h,:)');
        V_xx{n_h} = cost.phi_xx(x{n_h},x_ref(n_h,:)');
        
        % perform backwards propogation
        for k=n_h-1:-1:1
            [f_x, f_u] = model.f_d(x{k}, u{k}, dt);
            Q_x{k} = cost.L_x(x{k}, u{k}, x_ref(k,:)', u_ref(k,:)', dt) + ...
                     f_x.' * V_x{k+1};
                 
            Q_u{k} = cost.L_u(x{k}, u{k}, x_ref(k,:)', u_ref(k,:)', dt) + ...
                     f_u.' * V_x{k+1};
                 
            Q_xx{k} = cost.L_xx(x{k}, u{k}, x_ref(k,:)', u_ref(k,:)', dt) ...
                      + f_x.' * V_xx{k+1} * f_x;
                  
            Q_uu{k} = cost.L_uu(x{k}, u{k}, x_ref(k,:)', u_ref(k,:)', dt) ...
                      + f_u.' * V_xx{k+1} * f_u;
                  
            Q_xu{k} = cost.L_xu(x{k}, u{k}, x_ref(k,:)', u_ref(k,:)', dt) ...
                      + f_x.' * V_xx{k+1} * f_u;
                  
            Q_ux{k} = cost.L_ux(x{k}, u{k}, x_ref(k,:)', u_ref(k,:)', dt) ...
                      + f_u.' * V_xx{k+1} * f_x;
                  
            % Compute the value function derivatives
            V_x{k} = Q_x{k} - Q_xu{k} * (Q_uu{k} \ Q_u{k});
            V_xx{k} = Q_xx{k} - Q_xu{k} * (Q_uu{k} \ Q_ux{k});
        end
        
        % when to stop
        if i > 1
            if (J(i-1)-J(i)) <= stop_criterion 
                fprintf('[INFO]: Cost Converge at %d-th iteration with Cost of %7f.\n',i,J(i));
                break
            end
        end
    end
    SOL = u;
    K_f = K;
    x_rep = x;
end