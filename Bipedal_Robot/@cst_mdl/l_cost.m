function l = l_cost(in1,in2,in3,in4)
%L_COST
%    L = L_COST(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    19-May-2021 11:44:03

u1 = in2(1,:);
u2 = in2(2,:);
uref1 = in4(1,:);
uref2 = in4(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
xref1 = in3(1,:);
xref2 = in3(2,:);
xref3 = in3(3,:);
xref4 = in3(4,:);
xref5 = in3(5,:);
xref6 = in3(6,:);
et1 = ((u1./2.0-uref1./2.0).*(u1-uref1))./1.0e+3+((u2./2.0-uref2./2.0).*(u2-uref2))./1.0e+3+((x1./2.0e+3-xref1./2.0e+3).*(x1-xref1))./1.0e+3+((x2./2.0e+3-xref2./2.0e+3).*(x2-xref2))./1.0e+3+((x3./2.0e+3-xref3./2.0e+3).*(x3-xref3))./1.0e+3+((x4./2.0e+3-xref4./2.0e+3).*(x4-xref4))./1.0e+3;
et2 = ((x5./2.0e+3-xref5./2.0e+3).*(x5-xref5))./1.0e+3+((x6./2.0e+3-xref6./2.0e+3).*(x6-xref6))./1.0e+3;
l = et1+et2;
