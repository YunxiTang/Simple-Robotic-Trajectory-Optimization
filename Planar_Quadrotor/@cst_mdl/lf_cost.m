function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    30-Mar-2021 13:22:16

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
lf = (x2.*2.5e+1-2.5e+1).*(x2-1.0)+(x1.*2.5e+1-4.25e+2./2.0).*(x1-1.7e+1./2.0)+x3.^2.*2.5e+1+x4.^2.*2.5e+1+x5.^2.*2.5e+1+x6.^2.*2.5e+1;
