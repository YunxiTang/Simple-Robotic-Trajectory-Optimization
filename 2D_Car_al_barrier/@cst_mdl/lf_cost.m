function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    27-May-2021 18:08:55

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
lf = (x2.*(5.0./2.0)-1.5e+1./2.0).*(x2-3.0)+(x1.*(5.0./2.0)-2.5e+1./4.0).*(x1-5.0./2.0)+(x3-pi./2.0).*(x3.*(5.0./2.0)-pi.*(5.0./4.0))+x4.^2.*(5.0./2.0);
