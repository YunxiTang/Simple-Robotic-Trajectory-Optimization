function [V,V_x,V_xx,Q_x,Q_u,Q_uu,Q_xx,Q_xu,Q_ux,F_x,F_u,kff,Kfb] = ...
          backward_multi_prop(model, cost, Xtraj, Utraj, Xref, Uref, x_dft, dt, regType)
%BACKWARD_MULTI_PROP Do multiple shooting backward propogation
%%%        ARGS                                         DIMENSIONS
%%% INPUT: 
%%%        Xtraj            state traj(last iter)       [S, L+1, Nx, 1]
%%%        Utraj            control traj(last iter)     [S, L  , Nu, 1]
%%%        Xref             ref state traj(last iter)   [S, L+1, Nx, 1]
%%%        Uref             ref control traj(last iter) [S, L  , Nu, 1]
%%%        x_dft            defects along traj          [S, Nx, 1]
%%%        model:           robot model                 [class]
%%%        cost:            cost model                  [class]
%%%        umax:            max torque                  [Nu, 1]
%%%        dt:              grid time interval          [1]
%%%        regType:         Reguration Type             [1]
[S, LL , Nx, ~] = size(Xtraj);
[~, ~ , Nu, ~] = size(Utraj);
L = LL - 1;
N = S*L+1;
lambda = 1.0;
% Value function and derivatives
V = zeros(S, L+1, 1); 
V_x = zeros(S, L+1, Nx, 1); 
V_xx = zeros(S, L+1, Nx, Nx);

% State action value derivatives
Q_x = zeros(S, L, Nx, 1);   Q_u = zeros(S, L, Nu, 1);
Q_xx = zeros(S, L, Nx, Nx); Q_uu = zeros(S, L, Nu, Nu);
Q_xu = zeros(S, L, Nx, Nu); Q_ux = zeros(S, L, Nu, Nx);

F_x = zeros(S, L, Nx, Nx);
F_u = zeros(S, L, Nx, Nu);
kff = zeros(S, L, Nu, 1);
Kfb = zeros(S, L, Nx, 1);

% compute F_x, F_u
for i=1:S
    for j=1:L
        [f_x, f_u] = model.f_d(Xtraj(i,j,:,:), Utraj(i,j,:,:), dt);
        F_x(i,j,:,:) = f_x;
        F_u(i,j,:,:) = f_u;
    end
end

% compute terminal value function and derivatives
V(S,end,:) = cost.phi(Xtraj(S,L+1,:,:),Xref);
V_x(S,end,:,:) = cost.phi_x(Xtraj(S,L+1,:,:),Xref);
V_xx(S,end,:,:) = cost.phi_xx(Xtraj(S,L+1,:,:),Xref);

for i=S:-1:1
    for j=L:-1:1
        if j == L && 1 < i
            defect = x_dft(i-1,:,:)';
        else
            defect = zeros(Nx, 1);
        end
        fx = squeeze(F_x(i,j,:,:));
        fu = squeeze(F_u(i,j,:,:));
        
        Q_x(i,j,:,:) = cost.L_x(Xtraj(i,j,:,:), Utraj(i,j,:,:), Xref, zeros(Nu,1), dt) + ...
                                fx.' * (squeeze(V_x(i,j+1,:,:)) + squeeze(V_xx(i,j+1,:,:)) * defect);
                            
        Q_u(i,j,:,:) = cost.L_u(Xtraj(i,j,:,:), Utraj(i,j,:,:), Xref, zeros(Nu,1), dt) + ...
                                fu.' * (squeeze(V_x(i,j+1,:,:)) + squeeze(V_xx(i,j+1,:,:)) * defect);
        
        Q_xx(i,j,:,:) = cost.L_xx(Xtraj(i,j,:,:), Utraj(i,j,:,:), Xref, zeros(Nu,1), dt) + ...
                                  fx.' * squeeze(V_xx(i,j+1,:,:)) * fx;
                            
        Q_uu(i,j,:,:) = cost.L_uu(Xtraj(i,j,:,:), Utraj(i,j,:,:), Xref, zeros(Nu,1), dt) + ...
                                  fu.' * squeeze(V_xx(i,j+1,:,:)) * fu;
                              
        Q_ux(i,j,:,:) = cost.L_ux(Xtraj(i,j,:,:), Utraj(i,j,:,:), Xref, zeros(Nu,1), dt) + ...
                                  fu.' * squeeze(V_xx(i,j+1,:,:)) * fx;
                              
        Q_xu(i,j,:,:) = cost.L_xu(Xtraj(i,j,:,:), Utraj(i,j,:,:), Xref, zeros(Nu,1), dt) + ...
                                  fx.' * squeeze(V_xx(i,j+1,:,:)) * fu;
        
        V_xx_reg = squeeze(V_xx(i,j+1,:,:)) + lambda*eye(Nx)*(regType == 2);                      
                              
        Q_ux_reg = cost.L_ux(Xtraj(i,j,:,:), Utraj(i,j,:,:), Xref, zeros(Nu,1), dt) + ...
                             fu.' * V_xx_reg * fx;
        
        Q_uu_reg = cost.L_uu(Xtraj(i,j,:,:), Utraj(i,j,:,:), Xref, zeros(Nu,1), dt) + ...
                             fu.' * V_xx_reg * fu + lambda*eye(Nu)*(regType == 1);
                    
        l_n = -inv(Q_uu_reg) * squeeze(Q_u(i,j,:,:));
        L_n = -inv(Q_uu_reg) * Q_ux_reg;
        
        kff(i,j,:,:) = l_n;
        Kfb(i,j,:,:) = L_n;
        
        % Compute the value function derivatives
        V_x(i,j,:,:) = squeeze(Q_x(i,j,:,:)) - L_n' * Q_uu(i,j,:,:) * l_n;
        V_xx(i,j,:,:) = squeeze(Q_xx(i,j,:,:)) - L_n' * Q_uu(i,j,:,:) * L_n;
    end
    if i < S
        V(i,end,:) = V(i+1,1,:);
        V_x(i,end,:,:) = V_x(i+1,1,:,:);
        V_xx(i,end,:,:) = V_xx(i+1,1,:,:);
    end
end

end

