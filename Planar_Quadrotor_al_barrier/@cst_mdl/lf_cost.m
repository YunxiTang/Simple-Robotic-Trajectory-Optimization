function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    17-May-2021 18:50:36

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
lf = (x1.*2.5e+1-2.5e+1).*(x1-1.0)+(x2.*2.5e+1-7.5e+1./2.0).*(x2-3.0./2.0)+x3.^2.*2.5e+1+x4.^2.*2.5e+1+x5.^2.*2.5e+1+x6.^2.*2.5e+1;
