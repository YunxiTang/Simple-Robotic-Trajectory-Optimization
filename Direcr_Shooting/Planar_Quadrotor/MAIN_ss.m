clc; clear;
import casadi.*

params.T = 6;
params.N = 600;
params.x0 = [5.0; 2.5; 0.2; 0.0; 0.0; 0.0];
params.xf = [1.0; 1.5; 0.0; 0.0; 0.0; 0.0];
params.Q =  diag([1 1 1 1 1 1]);
params.R =  diag([0.1 0.1]);
params.Qf = diag([5 5 5 5 5 5])*10;
params.umax =  [5.0;5.0];
params.umin =  [0.1;0.1];
params.xub  = [];
params.xlb  = [];
params.nx = numel(params.x0);
params.nu = numel(params.umax);
% model variables
x = MX.sym('x',params.nx);
u = MX.sym('u',params.nu);

% model equations
xdot = Dynamics(x,u);

% objective function terms
L_path = 0.5*(x - params.xf)' * params.Q * (x - params.xf) + ...
         0.5*(u)' * params.R * (u);
L_final = 0.5*(x - params.xf)' * params.Qf * (x - params.xf);

% formulate the discrete dynamics
f = Function('f', {x, u}, {xdot, L_path});
qf = Function('qf', {x}, {L_final});

if false
    % CVODES from the SUNDIALS suite
%    dae = struct('x',x,'p',u,'ode',xdot,'quad',L_path);
%    opts = struct('tf',params.T/params.N);
%    F = integrator('F', 'cvodes', dae, opts);
else
    % Fixed step Runge-Kutta 4 integrator
    M = 1; % RK4 steps per interval
    DT = params.T/params.N/M;
    Q = 0;
    
    X0 = MX.sym('X0', params.nx);
    U = MX.sym('U',params.nu);
    X = X0;
    for j=1:M
       [k1, k1_q] = f(X, U);
       [k2, k2_q] = f(X + DT/2 * k1, U);
       [k3, k3_q] = f(X + DT/2 * k2, U);
       [k4, k4_q] = f(X + DT * k3, U);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
    end
    F = Function('F', {X0, U}, {X, Q}, {'x0','p'}, {'xf', 'qf'});
end
% Start with an empty NLP
w   = {};
w0  = [];
lbw = [];
ubw = [];
J   = 0;
g   = {};
lbg = [];
ubg = [];
%%
% Formulate the NLP
Xk = params.x0;
for k=0:(params.N-1)
    % New NLP variable for the control
    Uk = MX.sym(['U' num2str(k)],params.nu);
    w = {w{:}, Uk};
    lbw = [lbw; 
           params.umin];
    ubw = [ubw; 
           params.umax];
       
    w0 = [w0;  
         [3.924;3.924]];
    
    % Integrate till the end of the interval
    Fk = F('x0', Xk, 'p', Uk);
    Xk = Fk.xf;
    if k==(params.N-1)
        Jk = qf(Xk);
    else
        Jk = Fk.qf;
    end
    J = J + Jk;
    
    %%% Add inequality constraints
    gk = [0.4*0.4 - ((Xk(1)-2.0)*(Xk(1)-2.0)+(Xk(2)-2.0)*(Xk(2)-2.0));
          0.5*0.5 - ((Xk(1)-3.0)*(Xk(1)-3.0)+(Xk(2)-1.0)*(Xk(2)-1.0));
          0.4*0.4 - ((Xk(1)-4.0)*(Xk(1)-4.0)+(Xk(2)-2.0)*(Xk(2)-2.0));
          -Xk(2)-0;
           Xk(3)-deg2rad(30);
          -Xk(3)-deg2rad(30)];
    g = {g{:}, gk};
    lbg = [lbg; 
           [-inf;-inf;-inf;-inf;-inf;-inf]];
    ubg = [ubg; 
           [0;0;0;0;0;0]];

%     g = {g{:}, Xk};
%     lbg = [lbg; 
%            params.xlb];
%     ubg = [ubg; 
%            params.xub];
end

prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, ...
             'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

%% Plot the solution
u_opt = reshape(w_opt,[params.nu,params.N]);
figure();
plot(u_opt');
x_opt = params.x0;
Jr = 0.0;

for k=0:(params.N-1)
    Fk = F('x0',x_opt(:,end), 'p', u_opt(:,k+1));
    if k==(params.N-1)
        Jrk = qf(full(Fk.xf));
    else
        Jrk = Fk.qf;
    end
    Jr = Jr + full(Jrk);
    x_opt = [x_opt, full(Fk.xf)];
end
figure(2000);
xsol = x_opt;

figure();
plot(xsol');
legend('x1','x2','x3','x4');

%%
Obstacles = [2.0 3.0 4.0;
             2.0 1.0 2.0;
             0.4 0.5 0.4];
figure(2000);
plot_obstacle(Obstacles, 2000);hold on;
plot(params.x0(1), params.x0(2), 'kp', 'MarkerFaceColor', 'b', 'MarkerSize', 15); hold on;
plot(params.xf(1), params.xf(2), 'rh', 'MarkerFaceColor', 'r', 'MarkerSize', 15); hold on;
plot(xsol(1,:),xsol(2,:),'r-.','LineWidth',2.0);
h=legend('$Obstacle\;1$','$Obstacle\;2$','$Obstacle\;3$','$Start \; Point$','$Goal\;Point$','$Trajectory$', 'Interpreter','latex','FontSize',13);
h.NumColumns = 2;
set (gcf,'Position',[400,100,500,500], 'color','w');