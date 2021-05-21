function [success] = Setup_Functions_barrier(rbtmdl,path_constraintFunc,final_constraintFunc,obj_path_func, obj_bnd_func,params)
%SETUP_FUNCTIONS setup some important functions for solvers
%%% ARGS
%%% params           struct           set of params
%%% rbtmdl           class            robot model
%%% cstmdl           class            cost model
%%% success          int              1/0 (0->ERROR)
success = 1;
nx = params.nx;                         % dimension of state
nu = params.nu;                         % dimension of input
dt = sym('dt',[nu 1]','real');
u = sym('u',[nu 1]','real');            % symbolic variable: control
x = sym('x',[nx 1]','real');            % symbolic variable: state
xref = sym('xref',[nx 1]','real');      % symbolic variable: reference state
uref = sym('uref',[nu 1]','real');      % symbolic variable: reference input

%%% Relax-log barrier function
% path constraint
c = path_constraintFunc(x,u);
matlabFunction(c,  'vars', {x, u}, 'file', '@path_constraint/c',  'optimize',1==1);

cx = jacobian(c,x);
cu = jacobian(c,u);
matlabFunction(cx, 'vars', {x, u}, 'file', '@path_constraint/cx', 'optimize',1==1);
matlabFunction(cu, 'vars', {x, u}, 'file', '@path_constraint/cu', 'optimize',1==1);

% final state constraint
final_c = final_constraintFunc(x);
final_cx = jacobian(final_c, x);
matlabFunction(final_c,  'vars', {x}, 'file', '@final_constraint/final_c',  'optimize',1==1);
matlabFunction(final_cx, 'vars', {x}, 'file', '@final_constraint/final_cx', 'optimize',1==1);

% stage cost
l = obj_path_func(x, u, xref, uref);
lf = obj_bnd_func(x);

matlabFunction(l,  'vars', {x, u, xref, uref},...
               'file', '@cst_mdl/l_cost',  'optimize',1==1);

matlabFunction(lf, 'vars', {x}, 'file', '@cst_mdl/lf_cost', 'optimize',1==1);

lfx  = jacobian(lf, x)'; lfxx = jacobian(lfx,x);
matlabFunction(lf, lfx, lfxx, 'vars', {x}, 'file',...
              '@cst_mdl/lf_info', 'optimize',1==1);

%%% Note that jacobian(l,x) is actually a row vector when l is a scalar
lx   = jacobian(l,x)';  lu  = jacobian(l,u)';
lxx  = jacobian(lx,x);  luu  = jacobian(lu,u);
lux  = jacobian(lu,x);  lxu  = jacobian(lx,u);
matlabFunction(l,  lx,  lu, lxx, lux, lxu, luu,'vars',...
              {x, u, xref, uref}, 'file', '@cst_mdl/l_info','optimize', 1==1);
     
%%% derive  dynamics and compute anylytical derivatives
f = rbtmdl.Dynamics(0.0, x, u);
fx_continuous = jacobian(f,x);
fu_continuous = jacobian(f,u);
fx = eye(nx) + fx_continuous .* dt;
fu = fu_continuous .* dt;
matlabFunction(fx, fu,'vars',{x, u, dt},'file','@rbt_mdl/AnaLinSys', 'optimize', 1==1);

end

