function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    29-Apr-2021 10:32:34

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
lf = (x2.*1.25e+2-3.75e+2).*(x2-3.0)+(x1.*1.25e+2-6.25e+2./2.0).*(x1-5.0./2.0)+(x3-pi./2.0).*(x3.*1.25e+2-pi.*(1.25e+2./2.0))+x4.^2.*1.25e+2;
