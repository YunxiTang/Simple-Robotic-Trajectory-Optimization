 function [success] = Setup_Functions(params,rbtmdl,cstmdl,constraints)
%SETUP_FUNCTIONS setup some important functions for solvers
%%% ARGS
%%% params           struct           set of params
%%% rbtmdl           class            robot model
%%% cstmdl           class            cost model
%%% constraint       class            constraint structure
%%% success          int              1/0 (0->ERROR)
success = 0;
nx = params.nx;
nu = params.nu;
u = sym('u',[nu 1]','real');
x = sym('x',[nx 1]','real');
xref = sym('xref',[nx 1]','real');      % symbolic variable: reference state
uref = sym('uref',[nu 1]','real');      % symbolic variable: reference input

dt = params.dt;
Q = cstmdl.Q;
R = cstmdl.R;
Qf = cstmdl.Qf;

%%% write cost function (any form) here
%%% introduce relax-log barrier function later
l = 1/2*(x-xref).'*Q*(x-xref) + 1/2*(u-uref).'*R*(u-uref);
l = l * dt;
lf = 1/2*(x-xref).'*Qf*(x-xref);


%%% compute derivatives of cost function
%%% Note that jacobian(l,x) is actually a row vector
lx   = jacobian(l,x)';   lu  = jacobian(l,u)';
lxx  = jacobian(lx,x);   luu  = jacobian(lu,u);
lux  = jacobian(lu,x);   lxu  = jacobian(lx,u);
lfx  = jacobian(lf, x)'; lfxx = jacobian(lfx,x)';

matlabFunction(l,  'vars', {x, u, xref, uref}, 'file', '@cst_mdl/l_cost',  'optimize',1==1);
matlabFunction(lf, 'vars', {x, xref},          'file', '@cst_mdl/lf_cost', 'optimize',1==1);
matlabFunction(lf, lfx, lfxx, 'vars', {x, xref},'file', '@cst_mdl/lf_info', 'optimize',1==1);
matlabFunction(l,  lx,  lu, lxx, lux, lxu, luu, 'vars', {x, u, xref, uref}, 'file', '@cst_mdl/l_info','optimize', 1==1);

%%% derive dynamics and compute anylytical derivatives
f = rbtmdl.rk(x, u, dt);
fx = jacobian(f,x);
fu = jacobian(f,u);
matlabFunction(fx, fu,'vars',{x, u},'file','@planar_quadrotor/getLinSys', 'optimize',1==1);
success = 1;

%%% derive gradient of AL term from constraint 
c = constraints.c_ineq(x, u);
cx = jacobian(c, x);
cu = jacobian(c, u);
matlabFunction(cx, cu,'vars',{x, u},'file','@constraint/algrad', 'optimize',1==1);
end

