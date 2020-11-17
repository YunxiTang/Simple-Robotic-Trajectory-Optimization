%%% Local Stability Analysis for Van-der-Pol Oscillator
clear all;clc;yalmip('clear');
internal = zeros(2000,1);
%% General Settings
Iter_max = 5;
thershold = 1e-3;
BETA = [];
GAMMA = [];

sdpvar x1 x2 gamma
epsilon = 1e-5;
x = [x1;x2];
l = epsilon*x'*x;
ops = sdpsettings('solver','bmibnb','verbose',1,'debug',0,...
                  'bmibnb.lpsolver','linprog',...
                  'bmibnb.lowersolver','sedumi');
              
%% dynamics
f = [-x2;
      x1-x2+x2^3];
  
%% initial Search
% find the initial V0
S =[1.5000   -0.5000
   -0.5000    1.0000];

V0 = x'*S*x;
%% find gamma_star for V0
[s1,s1_c] = polynomial(x,2);
cons1 = -(l+jacobian(V0,x)*f)-s1*(gamma-V0);
prog = [sos(cons1),sos(s1),sos(gamma)];
solvesos(prog,-gamma,ops,[s1_c;gamma]);

disp('~~~~~~~~Gamma_star under the initial Lyapunov function is~~~~~~~~~~');
gamma_star = value(gamma);
disp(gamma_star);
% vanderplotter;

%% step 2: fix gamma_star and maximize the beta
ops1 = sdpsettings('solver','bmibnb','verbose',1,'debug',0,...
                  'bmibnb.lpsolver','linprog','bmibnb.lowersolver','sedumi',...
                  'bmibnb.uppersolver','fmincon',...
                  'bmibnb.maxiter',3000);
disp(gamma_star);
sdpvar beta

p = x'*[1 0;0 2]*x;
[s2,s2_c] = polynomial(x,2);
cons2 = (gamma_star-V0)-s2*(beta-p);
prog2 = [sos(cons2),sos(s2),sos(beta-epsilon)];
sol2 = solvesos(prog2,-beta,ops1,[s2_c;beta]);
disp('~~~~~~~~beta under the initial Lyapunov function is~~~~~~~~~~')
last_beta = value(beta);
disp(last_beta);

%% step 3: maximize beta by searching next higher order V (s1(x) and s2(x) fixed)
char_s1 = sdisplay(s1);char_s11 = char_s1{1}; % should be a string
char_s2 = sdisplay(s1);char_s22 = char_s2{1}; % should be a string
s1_c = value(s1_c); s2_c = value(s2_c);
yalmip('clear');
sdpvar x1 x2  beta
x = [x1;x2];
l = epsilon*x'*x;
p = x'*[1 0;0 2]*x;
f = [-x2;
      x1-x2+x2^3];
s1 = eval(char_s11);
s2 = eval(char_s22);

[V,V_c] = polynomial(x,4);

cons3 = -(l+jacobian(V,x)*f)-s1*(gamma_star-V);
cons4 = (gamma_star-V)-s2*(beta-p);
prog3 = [sos(V-l),sos(beta),sos(cons3),sos(cons4)];
sol3 = solvesos(prog3,-beta,[],[V_c;beta]);
disp('New 4th order V is found!');
disp(value(V_c));
BETA = [BETA last_beta];
GAMMA = [GAMMA gamma_star];

%% loop to iterative optimiztion
for i=1:Iter_max
    %% Fix V and get gamma_star under such a Lyaponov function and get gamma_star
    char_V = sdisplay(V);char_V = char_V{1};
    V_c = value(V_c); yalmip('clear');clear s1_c s2_c gamma_star
    sdpvar x1 x2 gamma
    x = [x1;x2];
    l = epsilon*x'*x;
    f = [-x2;
         x1-x2+x2^3];
    V = eval(char_V);
    disp('V is: ');
    sdisplay(V);
    %%
    [s1,s1_c] = polynomial(x,2);
    cons1 = -(l+jacobian(V,x)*f)-s1*(gamma-V);
    prog = [sos(cons1),sos(s1),sos(gamma)];
    solvesos(prog,-gamma,ops1,[s1_c;gamma]);

    disp('~~~~~~~~Gamma_star under the 4th Lyapunov function is~~~~~~~~~~');
    gamma_star = value(gamma);
    disp(gamma_star);
    GAMMA = [GAMMA gamma_star];
    
    %% Fix gamma_star, maximize beta
    sdpvar beta
    p = x'*[1 0;0 2]*x;
    [s2,s2_c] = polynomial(x,2);
    cons2 = (gamma_star-V)-s2*(beta-p);
    prog2 = [sos(cons2),sos(s2),sos(beta-epsilon)];
    sol2 = solvesos(prog2,-beta,ops1,[s2_c;beta]);
    disp('~~~~~~~~Beta under the 4th Lyapunov function is~~~~~~~~~~');
    last_beta = value(beta);
    disp(last_beta);
    BETA = [BETA last_beta];
    
    %% get next V function
    char_s1 = sdisplay(s1);char_s11 = char_s1{1}; 
    char_s2 = sdisplay(s1);char_s22 = char_s2{1}; 
    s1_c = value(s1_c); s2_c = value(s2_c);
    yalmip('clear'); clear s1 s2
    sdpvar x1 x2  beta
    x = [x1;x2];
    l = epsilon*x'*x;
    p = x'*[1 0;0 2]*x;
    f = [-x2;
          x1-x2+x2^3];
    s1 = eval(char_s11);
    s2 = eval(char_s22);

    [V,V_c] = polynomial(x,4);

    cons3 = -(l+jacobian(V,x)*f)-s1*(gamma_star-V);
    cons4 = (gamma_star-V)-s2*(beta-p);
    prog3 = [sos(V-l),sos(beta),sos(cons3),sos(cons4)];
    sol3 = solvesos(prog3,-beta,[],[V_c;beta]);
    disp('New 4th order V is found!');
    disp(value(V_c));
    
    fprintf('--------------------NEXT %d ITERATION-------------------',i);
    pause(1);
end






