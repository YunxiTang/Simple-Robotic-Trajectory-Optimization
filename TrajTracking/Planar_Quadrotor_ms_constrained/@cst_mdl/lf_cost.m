function lf = lf_cost(in1,in2)
%LF_COST
%    LF = LF_COST(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    18-May-2021 15:30:22

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
xref1 = in2(1,:);
xref2 = in2(2,:);
xref3 = in2(3,:);
xref4 = in2(4,:);
xref5 = in2(5,:);
xref6 = in2(6,:);
lf = (x1.*(5.0./2.0)-xref1.*(5.0./2.0)).*(x1-xref1)+(x2.*(5.0./2.0)-xref2.*(5.0./2.0)).*(x2-xref2)+(x3.*(5.0./2.0)-xref3.*(5.0./2.0)).*(x3-xref3)+(x4.*(5.0./2.0)-xref4.*(5.0./2.0)).*(x4-xref4)+(x5.*(5.0./2.0)-xref5.*(5.0./2.0)).*(x5-xref5)+(x6.*(5.0./2.0)-xref6.*(5.0./2.0)).*(x6-xref6);
