function l = l_cost(in1,in2)
%L_COST
%    L = L_COST(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Mar-2021 18:50:50

u1 = in2(1,:);
u2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
x4 = in1(4,:);
l = ((x1./2.0e+1-1.0./8.0).*(x1-5.0./2.0))./1.0e+2+((x2./2.0e+1-7.0./4.0e+1).*(x2-7.0./2.0))./1.0e+2+u1.^2./2.0e+3+u2.^2./2.0e+3+x4.^2./2.0e+3;
