clc; clear;
import casadi.*

params.T = 3;
params.N = 30;
params.xf = [0.0;pi;0.0;0.0];
params.Q = eye(4);
params.R = 0.1 * eye(1);
params.Qf = 500*eye(4);
params.umax = 25;
params.umin = -25;
params.xub  = [ 0.8; 1.5*pi; inf; inf];
params.xlb  = [-0.8;-1.5*pi;-inf;-inf];
    
% model variables
x = MX.sym('x',4,1);
u = MX.sym('u',1,1);

% model equations
xdot = Dynamics(x,u);

% objective function terms
L_path = (x - params.xf)' * params.Q * (x - params.xf) + ...
         (u)' * params.R * (u);
L_final = (x - params.xf)' * params.Qf * (x - params.xf);

% formulate the discrete dynamics
f = Function('f', {x, u}, {xdot, L_path});
qf = Function('qf', {x}, {L_final});
if false
    % CVODES from the SUNDIALS suite
   dae = struct('x',x,'p',u,'ode',xdot,'quad',L_path);
   opts = struct('tf',params.T/params.N);
   F = integrator('F', 'cvodes', dae, opts);
else
    % Fixed step Runge-Kutta 4 integrator
    M = 1; % RK4 steps per interval
    DT = params.T/params.N/M;
    Q = 0;
    
    X0 = MX.sym('X0', 4);
    U = MX.sym('U');
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
Xk = [0.0; deg2rad(0); 0.0; 0.0];
for k=0:(params.N-1)
    % New NLP variable for the control
    Uk = MX.sym(['U' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw, params.umin];
    ubw = [ubw, params.umax];
    w0 = [w0,  0];
    
    % Integrate till the end of the interval
    Fk = F('x0', Xk, 'p', Uk);
    Xk = Fk.xf;
    if k==(params.N-1)
        Jk = qf(Xk);
    else
        Jk = Fk.qf;
    end
    J = J + Jk;
    
    % Add inequality constraint
    g = {g{:}, Xk(1:2)};
    lbg = [lbg; 
           params.xlb(1:2)];
    ubg = [ubg; 
           params.xub(1:2)];
end

prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, ...
             'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

%% Plot the solution
u_opt = w_opt;
figure();
plot(u_opt);
x_opt = [0.0; deg2rad(0); 0.0; 0.0];
for k=0:(params.N-1)
    Fk = F('x0',x_opt(:,end), 'p', u_opt(k+1));
    x_opt = [x_opt, full(Fk.xf)];
end
figure();
plot(x_opt');
legend('x1','x2','x3','x4');