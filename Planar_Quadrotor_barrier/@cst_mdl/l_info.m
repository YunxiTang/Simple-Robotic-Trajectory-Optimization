function [l,lx,lu,lxx,lux,lxu,luu] = l_info(in1,in2,in3,in4)
%L_INFO
%    [L,LX,LU,LXX,LUX,LXU,LUU] = L_INFO(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    28-May-2021 14:00:23

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
l = ((u1./2.0e+1-uref1./2.0e+1).*(u1-uref1))./1.0e+2+((u2./2.0e+1-uref2./2.0e+1).*(u2-uref2))./1.0e+2+((x1.*2.0-xref1.*2.0).*(x1-xref1))./1.0e+2+((x2.*2.0-xref2.*2.0).*(x2-xref2))./1.0e+2+((x3.*2.0-xref3.*2.0).*(x3-xref3))./1.0e+2+((x4.*2.0-xref4.*2.0).*(x4-xref4))./1.0e+2+((x5.*2.0-xref5.*2.0).*(x5-xref5))./1.0e+2+((x6.*2.0-xref6.*2.0).*(x6-xref6))./1.0e+2;
if nargout > 1
    lx = [x1./2.5e+1-xref1./2.5e+1;x2./2.5e+1-xref2./2.5e+1;x3./2.5e+1-xref3./2.5e+1;x4./2.5e+1-xref4./2.5e+1;x5./2.5e+1-xref5./2.5e+1;x6./2.5e+1-xref6./2.5e+1];
end
if nargout > 2
    lu = [u1./1.0e+3-uref1./1.0e+3;u2./1.0e+3-uref2./1.0e+3];
end
if nargout > 3
    lxx = reshape([1.0./2.5e+1,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.5e+1,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.5e+1,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.5e+1,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.5e+1,0.0,0.0,0.0,0.0,0.0,0.0,1.0./2.5e+1],[6,6]);
end
if nargout > 4
    lux = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[2,6]);
end
if nargout > 5
    lxu = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,2]);
end
if nargout > 6
    luu = reshape([1.0./1.0e+3,0.0,0.0,1.0./1.0e+3],[2,2]);
end
