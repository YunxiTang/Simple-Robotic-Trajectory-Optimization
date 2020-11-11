%%% YALMIP for SOS
%% create symbolic decision variable
clear all;
x = sdpvar(1,1);
y = sdpvar(1,1);

p = (1+x)^4 + (1-y)^2;
F = sos(p);
solvesos(F);
h = sosd(F);
h_final = clean(p - h' * h,1e-6);

%% calculate a lower bound on a polynomial. (Parameterized problems)
clear all;
sdpvar x y lower_bound
p = (1+x*y)^2 - x*y + (1-y)^2;
F = sos(p - lower_bound);
solvesos(F,-lower_bound,[],lower_bound);
disp(value(lower_bound))

%% Multiple SOS constraints
clear all;
sdpvar x y t
p1 = t*(1+x*y)^2-x*y+(1-y)^2;
p2 = (1-x*y)^2+x*y+t*(1+y)^2;
F = [sos(p1),sos(p2)];
solvesos(F,t);
value(t)

%% Constrained polynomial optimization
sdpvar x y lower
p = (1+x*y)^2-x*y+(1-y)^2;
g = [1-x;1+x;1-y;1+y];

[s1,c1] = polynomial([x y],2);
[s2,c2] = polynomial([x y],2);
[s3,c3] = polynomial([x y],2);
[s4,c4] = polynomial([x y],2);

F = [sos(p-lower-[s1 s2 s3 s4]*g), sos(s1), sos(s2), sos(s3), sos(s4)];
solvesos(F,-lower,[],[c1;c2;c3;c4;lower]);

