function [fx,fu] = getLinSys(in1,u1)
%GETLINSYS
%    [FX,FU] = GETLINSYS(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    17-May-2021 10:24:10

x2 = in1(2,:);
x4 = in1(4,:);
t2 = cos(x2);
t3 = sin(x2);
t4 = x4.^2;
t5 = t2.^2;
t6 = t3.^2;
t7 = t5-5.0;
t8 = 1.0./t7;
t9 = t8.^2;
fx = reshape([1.0,0.0,0.0,0.0,0.0,1.0,t8.*(t5.*9.81e+2-t6.*9.81e+2+t2.*t4.*5.0e+1).*(-1.0e-4)-(t2.*t3.*t9.*(u1.*2.0e+2+t2.*t3.*9.81e+2+t3.*t4.*5.0e+1))./5.0e+3,(t8.*(t2.*9.81e+2+t4.*t5.*1.0e+1-t4.*t6.*1.0e+1-t3.*u1.*4.0e+1))./1.0e+3+(t2.*t3.*t9.*(t3.*9.81e+2+t2.*u1.*4.0e+1+t2.*t3.*t4.*1.0e+1))./5.0e+2,1.0./1.0e+2,0.0,1.0,0.0,0.0,1.0./1.0e+2,t3.*t8.*x4.*(-1.0./1.0e+2),(t2.*t3.*t8.*x4)./5.0e+1+1.0],[4,4]);
if nargout > 1
    fu = [0.0;0.0;t8.*(-1.0./5.0e+1);(t2.*t8)./2.5e+1];
end
