function [fx,fu] = getLinSys(in1,u1,noth1)
%GETLINSYS
%    [FX,FU] = GETLINSYS(IN1,U1,NOTH1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    23-Apr-2021 16:32:25

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = cos(x1);
t3 = cos(x2);
t4 = sin(x1);
t5 = sin(x2);
t6 = x1+x2;
t7 = x3.^2;
t8 = x4.^2;
t13 = x3.*6.6e+1;
t9 = t3.^2;
t10 = t5.^2;
t11 = cos(t6);
t12 = sin(t6);
t14 = t3.*1.0e+2;
t15 = -t13;
t16 = t5.*x3.*6.6e+3;
t17 = t5.*x4.*6.6e+3;
t18 = t3.*x3.*x4.*6.6e+3;
t20 = t2.*9.7119e+4;
t21 = t4.*9.7119e+4;
t23 = t3.*t8.*3.3e+3;
t24 = t5.*t8.*3.3e+3;
t27 = t3.*t5.*x3.*1.0e+4;
t28 = t3.*t5.*x4.*1.0e+4;
t19 = t16.*x4;
t22 = t9.*2.5e+3;
t25 = -t21;
t26 = t11.*1.30473e+5;
t29 = t5.*t12.*4.905e+4;
t30 = t3.*t11.*4.905e+4;
t31 = t3.*t12.*4.905e+4;
t32 = -t29;
t33 = t22-4.389e+3;
t34 = 1.0./t33;
t35 = t34.^2;
mt1 = [1.0,0.0,(t34.*(t20-t30))./2.0e+2,t34.*(t20-t26-t30+t2.*t3.*1.4715e+5).*(-1.0./2.0e+2),0.0,1.0,t34.*(t18+t23+t30+t32+t3.*t7.*3.3e+3+t7.*t9.*5.0e+3-t7.*t10.*5.0e+3+t5.*u1.*1.0e+4-t5.*x4.*1.0e+2).*(-1.0./2.0e+2)-t3.*t5.*t35.*(t15+t19+t24+t25+t31-u1.*6.6e+3+x4.*6.6e+1+t5.*t7.*3.3e+3-t3.*u1.*1.0e+4+t14.*x4+t3.*t5.*t7.*5.0e+3).*2.5e+1];
mt2 = [(t34.*(t18+t23+t26+t30+t32+t4.*t5.*1.4715e+5+t3.*t7.*1.66e+4+t7.*t9.*1.0e+4-t7.*t10.*1.0e+4+t8.*t9.*5.0e+3-t8.*t10.*5.0e+3+t5.*u1.*2.0e+4+t5.*x3.*1.0e+2-t5.*x4.*2.0e+2+t9.*x3.*x4.*1.0e+4-t10.*x3.*x4.*1.0e+4))./2.0e+2+t3.*t5.*t35.*(t12.*1.30473e+5+t15+t19+t24+t25+t31-u1.*3.32e+4+x4.*3.32e+2-t3.*t4.*1.4715e+5+t5.*t7.*1.66e+4-t3.*u1.*2.0e+4-t3.*x3.*1.0e+2+t3.*x4.*2.0e+2+t27.*x4+t3.*t5.*t7.*1.0e+4+t3.*t5.*t8.*5.0e+3).*2.5e+1];
mt3 = [1.0./1.0e+2,0.0,t34.*(t16+t17+t27-6.6e+1).*(-1.0./2.0e+2)+1.0,(t34.*(-t14+t17+t28+t5.*x3.*3.32e+4+t3.*t5.*x3.*2.0e+4-6.6e+1))./2.0e+2,0.0,1.0./1.0e+2,t34.*(t14+t16+t17+6.6e+1).*(-1.0./2.0e+2),(t34.*(t3.*2.0e+2+t16+t17+t27+t28+3.32e+2))./2.0e+2+1.0];
fx = reshape([mt1,mt2,mt3],4,4);
if nargout > 1
    fu = [0.0;0.0;(t34.*(t3.*1.0e+4+6.6e+3))./2.0e+2;t34.*(t3.*2.0e+4+3.32e+4).*(-1.0./2.0e+2)];
end
