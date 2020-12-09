function [SOL] = DDP(x_0,x_star,t_f,N,dyn,cost,u_max,num_iter,alpha,stop_criterion)
%DDP DDP process
    %% allocate arrays memory for DDP structure
    % Time stamp array and timestep
    t = linspace(0.0,t_f,N);
    dt = t(2) - t(1);
    
    % cost history
    J = zeros(num_iter,1);
    % control energy history
    Eu = zeros(num_iter,1);
    
    % state history
    x = cell(N,1);
    x_new = cell(N,1);
    
    x{1} = x_0;
    x_new{1} = x_0;
    
    % control input his
    u = cell(N,1);
    
    % Value function and derivatives
    V = zeros(N,1);
    V_x = cell(N,1);
    V_xx = cell(N,1);
    
    % State action value derivatives
    Q_x = cell(N,1);
    Q_u = cell(N,1);
    Q_xx = cell(N,1);
    Q_uu = cell(N,1);
    Q_xu = cell(N,1);
    Q_ux = cell(N,1);
    
    %% initialize DDP with random input sequence
    disp('Initializing the DDP with random control sequence.....');
    
    % generate the random input
    u{N} = zeros(numel(u_max),1);
    for i=1:N-1
        u{i} = 2 * u_max .* rand(numel(u_max),1) - u_max;
    end
    
    % generate initial trajectory using random control input (Forward Simulation)
    disp('Generating the initial trajectory (RK45).');
    for j=1:N-1
        % RK45 Integration
        k1 = dyn.F(x_new{j},             u{j});
        k2 = dyn.F(x_new{j} + 0.5*dt*k1, u{j});
        k3 = dyn.F(x_new{j} + 0.5*dt*k2, u{j});
        k4 = dyn.F(x_new{j} +     dt*k3, u{j});
        x_new{j+1} = x_new{j} + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4);
    end

    %% perform main DDP iterations on the trajectory and input sequence
    disp('Begin DDP iterations...');
    
    % get into the loop of iteration (DDP).
    for i=1:num_iter
        fprintf(' -------DDP Iteration %d out of %d------- \n',i,num_iter);
        
        % FORWARD SIMULATION WITH updated control sequence from previous iteration
        if i > 1
            for k=1:N-1
                % compute control update feed-forward and feed-back
                du_ff = -inv(Q_uu{k}) * Q_u{k};
                du_fb = -inv(Q_uu{k}) * Q_ux{k} * (x_new{k} - x{k});
                
                % limit feed-forward control modification with clamping
                for m = 1:numel(u_max)
                    du_ff(m) = min(u_max(m), max(-u_max(m), ...
                                             du_ff(m) + u{k}(m))) - u{k}(m);
                end
            
                % update control
                u{k} = u{k} + alpha .* (du_ff) + du_fb;
                
                % FORWARD SIMULATION: Compute next state in trajectory with new control
                k1 = dyn.F(x_new{k},             u{k});
                k2 = dyn.F(x_new{k} + 0.5*dt*k1, u{k});
                k3 = dyn.F(x_new{k} + 0.5*dt*k2, u{k});
                k4 = dyn.F(x_new{k} +     dt*k3, u{k});
                x_new{k+1} = x_new{k} + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4);
                
                % Return error if problem with trajectory
                if isnan(x_new{k+1})
                    SOL = assemble_solution(x, u, t, J, Eu, Q_u, ...
                                            Q_uu, Q_ux, 1);
                    return
                end
            end 
        end
        
        % update trajectory
        x = x_new;
        
        % compute total cost
        J(i) = cost.phi(x{N},x_star);
        
        for k=1:N-1
            J(i) = J(i) + cost.L(x{k}, u{k}, dt);
        end
        
        % compute control effort
        for k = 1:N-1
            Eu(i) = Eu(i) + 0.5 * u{k}.' * u{k} * dt; 
        end
        
        % compute terminal value function and derivatives
        V(N) = cost.phi(x{N},x_star);
        V_x{N} = cost.phi_x(x{N},x_star);
        V_xx{N} = cost.phi_xx(x{N},x_star);
        
        % perform backwards propogation
        for k=N-1:-1:1
            Q_x{k} = cost.L_x(x{k}, u{k}, dt) + ...
                     dyn.Phi(x{k},u{k}, dt).' * V_x{k+1};
            Q_u{k} = cost.L_u(x{k}, u{k}, dt) + ...
                     dyn.beta(x{k}, u{k}, dt).' * V_x{k+1};
            Q_xx{k} = cost.L_xx(x{k}, u{k}, dt) ...
                      + dyn.Phi(x{k}, u{k}, dt).' * V_xx{k+1} ...
                      * dyn.Phi(x{k}, u{k}, dt);
            Q_uu{k} = cost.L_uu(x{k}, u{k}, dt) ...
                      + dyn.beta(x{k}, u{k}, dt).' * V_xx{k+1} ...
                      * dyn.beta(x{k}, u{k}, dt);
            Q_xu{k} = cost.L_xu(x{k}, u{k}, dt) ...
                      + dyn.Phi(x{k}, u{k}, dt).' * V_xx{k+1} ...
                      * dyn.beta(x{k}, u{k}, dt);
            Q_ux{k} = cost.L_ux(x{k}, u{k}, dt) ...
                      + dyn.beta(x{k}, u{k}, dt).' * V_xx{k+1} ...
                      * dyn.Phi(x{k}, u{k}, dt);
                  
            % Compute the value function derivatives
            V_x{k} = Q_x{k} - Q_xu{k} * (Q_uu{k} \ Q_u{k});
            V_xx{k} = Q_xx{k} - Q_xu{k} * (Q_uu{k} \ Q_ux{k});
        end
        
        % when to stop
        if i>1
            if abs(J(i)-J(i-1)) <= stop_criterion
                fprintf('~~~~~~Cost Converge at %d-th iteration with Cost of %f. ~~~~~~~~\n',i,J(i));
                break
            end
        end
    end
    %% Assemble and return solution structure
    if i ==  num_iter
        fprintf('Reach the Maximun Iteration Number %d.\n', num_iter);
    end
    fprintf("finished DDP, assembling results for post-processing...\n");
    
    % Assemble solution
    SOL = assemble_solution(x, u, t, J, Eu, Q_u, Q_uu, Q_ux, 0);
end


function sol = assemble_solution(x, u, t, J, E, Q_u, Q_uu, Q_ux, error)

    % Solution structure
    sol = struct;
    sol.error = error;
    sol.x = x;
    sol.u = u;
    sol.t = t;
    sol.dt = t(2) - t(1);
    sol.J = J;
    sol.E = E;
    sol.Q_u = Q_u;
    sol.Q_uu = Q_uu;
    sol.Q_ux = Q_ux;
    
    return 
end
