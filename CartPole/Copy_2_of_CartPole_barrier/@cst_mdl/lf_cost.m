function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    21-May-2021 16:18:53

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
lf = (x2-pi).*(x2.*(5.0./2.0)-pi.*(5.0./2.0))+x1.^2.*(5.0./2.0)+x3.^2.*(5.0./2.0)+x4.^2.*(5.0./2.0);
