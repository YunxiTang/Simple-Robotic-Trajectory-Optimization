clear all;
pvar x1 x2;
vartable = [x1;x2];
% initialize the SOSP
prog = sosprogram(vartable);
    
% inequality
% p(x1,x2)>=0
p = 2*x1^4 + 2*x1^3*x2 - x1^2*x2^2 + 5*x2^4;
prog = sosineq(prog,p);
solver_opt.solver = 'sedumi';
[SOSP,INFO] = sossolve(prog); % if feasible, then p is SOS


