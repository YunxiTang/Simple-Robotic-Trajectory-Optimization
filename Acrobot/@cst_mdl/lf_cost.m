function lf = lf_cost(in1)
%LF_COST
%    LF = LF_COST(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    28-Apr-2021 21:00:52

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
lf = (x1.*6.0e+1-9.42e+2./5.0).*(x1-1.57e+2./5.0e+1)+x2.^2.*6.0e+1+x3.^2.*6.0e+1+x4.^2.*6.0e+1;
