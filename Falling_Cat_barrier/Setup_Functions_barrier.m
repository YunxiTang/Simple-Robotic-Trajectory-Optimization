function [success] = Setup_Functions_barrier(params,rbtmdl,cstmdl,cons)
%SETUP_FUNCTIONS setup some important functions for solvers
%%% ARGS
%%% params           struct           set of params
%%% rbtmdl           class            robot model
%%% cstmdl           class            cost model
%%% success          int              1/0 (0->ERROR)
success = 0;
nx = params.nx;
nu = params.nu;
u = sym('u',[nu 1]','real');
x = sym('x',[nx 1]','real');
xref = sym('xref',[nx 1]','real');
uref = sym('uref',[nu 1]','real');

xf = params.xf;
dt = params.dt;
Q = cstmdl.Q;
R = cstmdl.R;
Qf = cstmdl.Qf;

%%% write cost function (any form) here
%%% introduce relax-log barrier function later
c = cons(x,u);
cx = jacobian(c,x);
cu = jacobian(c,u);
matlabFunction(c, 'vars',{x, u},'file','@constraint/c', 'optimize',1==1);
matlabFunction(cx, 'vars',{x, u},'file','@constraint/cx', 'optimize',1==1);
matlabFunction(cu, 'vars',{x, u},'file','@constraint/cu', 'optimize',1==1);

l = 1/2*(x-xref).'*Q*(x-xref) + 1/2*(u-uref).'*R*(u-uref);
l = l * dt;
lf = 1/2*(x-xf).'*Qf*(x-xf);


%%% compute derivatives of cost function
%%% Note that jacobian(l,x) is actually a row vector
lx   = jacobian(l,x)';  lu  = jacobian(l,u)';
lxx  = jacobian(lx,x);  luu  = jacobian(lu,u);
lux  = jacobian(lu,x);  lxu  = jacobian(lx,u);
lfx  = jacobian(lf, x); lfxx = jacobian(lfx,x);
matlabFunction(l, 'vars',{x, u, xref, uref},'file','@cst_mdl/l_cost', 'optimize',1==1);
matlabFunction(lf,'vars',{x},'file','@cst_mdl/lf_cost','optimize',1==1);
matlabFunction(l,lx,lu,lxx,lux,lxu,luu, 'vars',{x, u, xref, uref},...
               'file','@cst_mdl/l_info','optimize',1==1);
matlabFunction(lf,lfx,lfxx,'vars',{x},'file','@cst_mdl/lf_info','optimize',1==1);

%%% derive  dynamics and compute anylytical derivatives
% f = rbtmdl.Dynamics(0.0, x, u);
% fx_continuous = jacobian(f,x);
% fu_continuous = jacobian(f,u);
% fx = eye(nx) + fx_continuous .* dt;
% fu = fu_continuous .* dt;
% matlabFunction(fx, fu,'vars',{x, u},'file','@falling_cat/getLinSys', 'optimize',1==0);
end

