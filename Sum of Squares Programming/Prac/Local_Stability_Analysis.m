%%% Local Stability Analysis for Van-der-Pol Oscillator
clear all;clc;
pvar x y
mu=1; r=2.819835;
g = r - (x^2 + y^2);
f = [-y; 
     -mu * (1 - x^2 ) * y + x];
prog=sosprogram([x y]);
Z2=monomials([x y],1:2);
Z4=monomials([x y],0:4);
[prog,V]=sossosvar(prog,Z2,'wscoeff');
V = V + .0001 * (x^4 + y^4);

nablaV=[diff(V,x);diff(V,y)];
[prog,s]=sossosvar(prog,Z2);
prog=sosineq(prog,-nablaV'*f-s*g);
prog=sossolve(prog);
Vn=sosgetsol(prog,V);