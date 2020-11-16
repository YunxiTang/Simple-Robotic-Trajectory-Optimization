%%% Local Stability Analysis for Van-der-Pol Oscillator
clear all;clc;yalmip('clear')
sdpvar x1 x2 gamma
epsilon = 1e-5;
x = [x1;x2];
l = epsilon*x'*x;
ops = sdpsettings('solver','bmibnb','verbose',1,'debug',0,...
                  'bmibnb.lpsolver','linprog','bmibnb.lowersolver','sedumi');
%% dynamics
f = [-x2;
      x1-x2+x2^3];

%% find the initial V0
S =[1.5000   -0.5000
   -0.5000    1.0000];

V0 = x'*S*x;

%%
[s1,s1_c] = polynomial(x,2);

cons1 = -(l+jacobian(V0,x)*f)-s1*(gamma-V0);

prog = [sos(cons1),sos(s1),sos(gamma)];

solvesos(prog,-gamma,ops,[s1_c;gamma]);

%%
disp('~~~~~~~~Gamma under the initial Lyapunov function is~~~~~~~~~~');
value(gamma);
gamma_star = value(gamma);
vanderplotter;
%% step 2
ops1 = sdpsettings('solver','bmibnb','verbose',1,'debug',1,...
                  'bmibnb.lpsolver','linprog','bmibnb.lowersolver','sedumi',...
                  'bmibnb.uppersolver','fmincon',...
                  'bmibnb.maxiter',1500);
disp(gamma_star);
sdpvar beta

p = x'*x/2;
[s2,s2_c] = polynomial(x,2);
cons2 = (gamma_star-V0)-s2*(beta-p);
prog2 = [sos(cons2),sos(s2),sos(beta)];
sol2 = solvesos(prog2,-beta,ops1,[s2_c;beta]);
ApproxROA;

%% step 3
yalmip('clear')
sdpvar x1 x2 beta
epsilon = 1e-5;
x = [x1;x2];
l = epsilon*x'*x;

t1 = sdisplay(s1);
t2 = char(t1);
s1_c = value(s1_c);
s1 = eval(t2);


[V,V_c] = polynomial(x,4);
cons3 = -(l+jacobian(V,x)*f)-s1*(gamma_star-V);
cons4 = (gamma_star-V)-s2*(beta-p);
prog3 = [sos(V-l),sos(beta),sos(cons3),sos(cons4)];
sol3 = solvesos(prog3,-beta,[],[V_c;beta]);


