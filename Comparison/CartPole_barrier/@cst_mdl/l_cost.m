function l = l_cost(in1,u1,in3,uref1)
%L_COST
%    L = L_COST(IN1,U1,IN3,UREF1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    06-May-2021 11:18:23

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
xref1 = in3(1,:);
xref2 = in3(2,:);
xref3 = in3(3,:);
xref4 = in3(4,:);
l = ((u1./2.0e+1-uref1./2.0e+1).*(u1-uref1))./1.0e+2+((x1./2.0-xref1./2.0).*(x1-xref1))./1.0e+2+((x2./2.0-xref2./2.0).*(x2-xref2))./1.0e+2+((x3./2.0-xref3./2.0).*(x3-xref3))./1.0e+2+((x4./2.0-xref4./2.0).*(x4-xref4))./1.0e+2;
