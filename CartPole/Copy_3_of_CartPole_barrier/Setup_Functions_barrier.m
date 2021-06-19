function [success] = Setup_Functions_barrier(params,rbtmdl,cstmdl,path_constraint,final_constraint)
%SETUP_FUNCTIONS setup some important functions for solvers
%%% ARGS
%%% params           struct           set of params
%%% rbtmdl           class            robot model
%%% cstmdl           class            cost model
%%% success          int              1/0 (0->ERROR)
success = 1;
nx = params.nx;                         % dimension of state
nu = params.nu;                         % dimension of input
u = sym('u',[nu 1]','real');            % symbolic variable: control
x = sym('x',[nx 1]','real');            % symbolic variable: state
xref = sym('xref',[nx 1]','real');      % symbolic variable: reference state
uref = sym('uref',[nu 1]','real');      % symbolic variable: reference input
dt   = sym('dt','real');  
xf = params.xf;
% dt = params.dt;
Q = cstmdl.Q;
R = cstmdl.R;
Qf = cstmdl.Qf;

%%% Relax-log barrier function
% path constraint
c = path_constraint(x,u);
cx = jacobian(c,x);
cu = jacobian(c,u);
matlabFunction(c,  'vars', {x, u}, 'file', '@path_constraint/c',  'optimize',1==1);
matlabFunction(cx, 'vars', {x, u}, 'file', '@path_constraint/cx', 'optimize',1==1);
matlabFunction(cu, 'vars', {x, u}, 'file', '@path_constraint/cu', 'optimize',1==1);

% final state constraint
final_c = final_constraint(x);
final_cx = jacobian(final_c, x);
matlabFunction(final_c,  'vars', {x}, 'file', '@final_constraint/final_c',  'optimize',1==1);
matlabFunction(final_cx, 'vars', {x}, 'file', '@final_constraint/final_cx', 'optimize',1==1);

% stage cost
l = 1/2*(x-xref).'*Q*(x-xref) + 1/2*(u-uref).'*R*(u-uref);
l = l * params.dt;
lf = 1/2*(x-xf).'*Qf*(x-xf);


%%% compute derivatives of cost function
%%% Note that jacobian(l,x) is actually a row vector when l is a scalar
lx   = jacobian(l,x)';  lu  = jacobian(l,u)';
lxx  = jacobian(lx,x);  luu  = jacobian(lu,u);
lux  = jacobian(lu,x);  lxu  = jacobian(lx,u);
lfx  = jacobian(lf, x)'; lfxx = jacobian(lfx,x);

matlabFunction(l,  'vars', {x, u, xref, uref}, 'file', '@cst_mdl/l_cost',  'optimize',1==1);
matlabFunction(lf, 'vars', {x},                'file', '@cst_mdl/lf_cost', 'optimize',1==1);
matlabFunction(lf, lfx, lfxx, 'vars', {x},     'file', '@cst_mdl/lf_info', 'optimize',1==1);
matlabFunction(l,  lx,  lu, lxx, lux, lxu, luu, 'vars', {x, u, xref, uref}, 'file', '@cst_mdl/l_info','optimize', 1==1);

%%% derive  dynamics and compute anylytical derivatives
f = rbtmdl.Dynamics(0.0, x, u);
fx_continuous = jacobian(f,x);
fu_continuous = jacobian(f,u);
fx = eye(nx) + fx_continuous .* dt;
fu = fu_continuous .* dt;
rbtmdl.getLinSys = matlabFunction(fx, fu,'vars',{x, u, dt}, 'optimize',1==0);
end

