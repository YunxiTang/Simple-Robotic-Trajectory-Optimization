function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    16-Apr-2021 14:59:56

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x6 = in1(6,:);
x7 = in1(7,:);
x8 = in1(8,:);
lf = x1.^2.*3.0+x2.^2.*3.0+x3.^2.*3.0+x6.^2.*3.0+x7.^2.*3.0+x8.^2.*3.0;
