function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    13-Apr-2021 14:40:49

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
lf = (x1.*2.5e+2-7.0e+2).*(x1-1.4e+1./5.0)+(x2.*2.5e+2-7.0e+2).*(x2-1.4e+1./5.0)+(x3-pi./2.0).*(x3.*2.5e+2-pi.*1.25e+2)+x4.^2.*2.5e+2;
