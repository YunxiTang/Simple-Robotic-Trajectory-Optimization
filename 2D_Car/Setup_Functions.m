function [success] = Setup_Functions(params,rbtmdl,cstmdl)
%SETUP_FUNCTIONS setup some important functions for solvers
%%% ARGS
%%% params           struct           set of params
%%% rbtmdl           class            robot model
%%% cstmdl           class            cost model
%%% success          int              1/0 (0->ERROR)
success = 0;
alpha = 1e-5;
nx = params.nx;
nu = params.nu;
u = sym('u',[nu 1]','real');
x = sym('x',[nx 1]','real');

xf = params.xf;
dt = params.dt;
Q = cstmdl.Q;
R = cstmdl.R;
Qf = cstmdl.Qf;

%%% write cost function (any form) here
%%% introduce relax-log barrier function later

l = sqrt((x).'*Q*(x)+alpha^2)-alpha + 0.5 * u.'*R*u;
l = l * dt;
lf = sqrt((x-xf).'*Qf*(x-xf)+alpha^2)-alpha;


%%% compute derivatives of cost function
%%% Note that jacobian(l,x) is actually a row vector
lx   = jacobian(l,x)';  lu  = jacobian(l,u)';
lxx  = jacobian(lx,x);  luu  = jacobian(lu,u);
lux  = jacobian(lu,x);  lxu  = jacobian(lx,u);
lfx  = jacobian(lf, x); lfxx = jacobian(lfx,x);
matlabFunction(l, 'vars',{x, u},'file','@cst_mdl/l_cost', 'optimize',1==1);
matlabFunction(lf,'vars',{x},'file','@cst_mdl/lf_cost','optimize',1==1);
matlabFunction(l,lx,lu,lxx,lux,lxu,luu, 'vars',{x, u},...
               'file','@cst_mdl/l_info','optimize',1==1);
matlabFunction(lf,lfx,lfxx,'vars',{x},'file','@cst_mdl/lf_info','optimize',1==1);

%%% derive Car dynamics and compute anylytical derivatives
f = rbtmdl.rk(x, u, dt);
fx = jacobian(f,x);
fu = jacobian(f,u);
matlabFunction(fx, fu,'vars',{x, u},'file','@Car/getLinSys', 'optimize',1==1);
success = 1;
end

