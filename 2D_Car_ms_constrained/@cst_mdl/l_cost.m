function l = l_cost(in1,in2)
%L_COST
%    L = L_COST(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    27-Apr-2021 15:56:53

u1 = in2(1,:);
u2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
l = ((x2./2.0-3.0./2.0).*(x2-3.0))./1.0e+2+((x1./2.0-5.0./4.0).*(x1-5.0./2.0))./1.0e+2+((x3-pi./2.0).*(x3./2.0-pi./4.0))./1.0e+2+u1.^2./2.0e+3+u2.^2./2.0e+3+x4.^2./2.0e+2;
