function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Mar-2021 21:56:10

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
lf = (x1.*1.0e+1-4.0e+1).*(x1-4.0)+(x2.*1.0e+1-4.0e+1).*(x2-4.0)+x3.^2.*1.0e+1+x4.^2.*1.0e+1;
