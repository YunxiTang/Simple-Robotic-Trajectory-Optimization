 function [success] = Setup_Functions(params,rbtmdl,cstmdl,path_constraints,final_constraints)
%SETUP_FUNCTIONS setup some important functions for solvers
%%% ARGS
%%% params           struct           set of params
%%% rbtmdl           class            robot model
%%% cstmdl           class            cost model
%%% constraint       class            constraint structure
%%% success          int              1/0 (0->ERROR)
success = 0;
alpha = 1e-1;
nx = params.nx;
nu = params.nu;
u = sym('u',[nu 1]','real');
x = sym('x',[nx 1]','real');
internal = sym('internal',[1 1]','real');

xf = params.xf;
dt = params.dt;
Q = cstmdl.Q;
R = cstmdl.R;
Qf = cstmdl.Qf;

%%% write cost function (any form) here
l = 1/2*(x-xf).'*Q*(x-xf) + 1/2*u.'*R*u;
% l = sqrt((x-xf).'*Q*(x-xf)+alpha^2)-alpha + 0.5 * u.'*R*u;
l = l * dt;
lf = 1/2*(x-xf).'*Qf*(x-xf);


%%% compute derivatives of cost function
%%% Note that jacobian(l,x) is actually a row vector
lx   = jacobian(l,x)';   lu  = jacobian(l,u)';
lxx  = jacobian(lx,x);   luu  = jacobian(lu,u);
lux  = jacobian(lu,x);   lxu  = jacobian(lx,u);
lfx  = jacobian(lf, x)'; lfxx = jacobian(lfx,x);

matlabFunction(l, 'file','@cst_mdl/l_cost', 'vars',{x, u}, 'outputs',{'l'});
matlabFunction(lf,'file','@cst_mdl/lf_cost','vars',{x},'outputs',{'lf'});
matlabFunction(l,lx,lu,lxx,lux,lxu,luu, 'file','@cst_mdl/l_info', 'vars',{x, u}, 'outputs', {'l','lx','lu','lxx','lux','lxu','luu'});
matlabFunction(lf,lfx,lfxx,'file','@cst_mdl/lf_info','vars',{x},'outputs', {'lf','lfx','lfxx'});
disp('[INFO]: Cost-related Function Generation done.');

%%% derive gradient of AL term from constraint 
c = path_constraints.c(x, u);
cx = jacobian(c, x);
cu = jacobian(c, u);
matlabFunction(cx, 'file','@path_constraint/cx', 'vars',{x, u},'outputs',{'cx'});
matlabFunction(cu, 'file','@path_constraint/cu', 'vars',{x, u},'outputs',{'cu'});
disp('[INFO]: Constraint-related Function Generation done.');

cf= final_constraints.c(x);
cfx = jacobian(cf, x);
matlabFunction(cfx, 'file','@final_constraint/cx', 'vars',{x, u},'outputs',{'cx'});

disp('[INFO]: Constraint-related Function Generation done.');

%%% derive dynamics and compute anylytical derivatives
% f = rbtmdl.Dynamics(0, x, u);
% fx = eye(nx) + jacobian(f,x) * internal;
% fu = jacobian(f,u) * internal;
% matlabFunction(fx, fu, 'file','@planararm/getLinSys', 'vars',{x, u, internal}, 'outputs', {'fx','fu'});
% disp('[INFO]: Dynamics-related Function Generation done.');
end

