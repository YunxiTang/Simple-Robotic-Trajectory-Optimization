function [fx,fu] = getLinSys(in1,in2)
%GETLINSYS
%    [FX,FU] = GETLINSYS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    20-May-2021 16:41:48

u1 = in2(1,:);
u2 = in2(2,:);
x3 = in1(3,:);
x6 = in1(6,:);
t2 = cos(x3);
t3 = sin(x3);
t4 = u1.*(5.0./4.0);
t5 = u2.*(5.0./4.0);
t6 = x6./1.0e+2;
t7 = x6./2.0e+2;
t10 = u1./1.92e+3;
t11 = u2./1.92e+3;
t12 = u1./3.84e+3;
t13 = u2./3.84e+3;
t8 = t2./4.8e+2;
t9 = t3./4.8e+2;
t15 = t7+x3;
t16 = -t10;
t17 = -t12;
t20 = t2./4.8e+4;
t21 = t3./4.8e+4;
t22 = t4+t5;
t14 = -t9;
t18 = cos(t15);
t19 = sin(t15);
t23 = -t21;
t30 = t6+t11+t16+x3;
t31 = t13+t15+t17;
t24 = t18./2.4e+2;
t25 = t19./2.4e+2;
t27 = t18./4.8e+4;
t28 = t19./4.8e+4;
t32 = cos(t30);
t33 = cos(t31);
t34 = sin(t30);
t35 = sin(t31);
t36 = (t18.*t22)./6.0e+4;
t37 = (t19.*t22)./6.0e+4;
t26 = -t25;
t29 = -t28;
t38 = -t36;
t39 = -t37;
t40 = t32./4.8e+2;
t41 = t33./2.4e+2;
t42 = t34./4.8e+2;
t43 = t35./2.4e+2;
t46 = t33./4.8e+4;
t47 = t35./4.8e+4;
t49 = (t22.*t33)./6.0e+4;
t50 = (t22.*t35)./6.0e+4;
t53 = (t22.*t32)./1.152e+6;
t54 = (t22.*t33)./1.152e+6;
t55 = (t22.*t34)./1.152e+6;
t56 = (t22.*t35)./1.152e+6;
t57 = (t22.*t33)./2.304e+8;
t58 = (t22.*t35)./2.304e+8;
t44 = -t42;
t45 = -t43;
t48 = -t47;
t51 = -t49;
t52 = -t50;
mt1 = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,t38+t51-(t2.*t22)./6.0e+4,t39+t52-(t3.*t22)./6.0e+4,1.0,t2.*t22.*(-1.0./6.0e+2)-(t18.*t22)./3.0e+2-(t22.*t32)./6.0e+2-(t22.*t33)./3.0e+2,t3.*t22.*(-1.0./6.0e+2)-(t19.*t22)./3.0e+2-(t22.*t34)./6.0e+2-(t22.*t35)./3.0e+2,0.0,1.0./1.0e+2,0.0,0.0,1.0,0.0,0.0,0.0,1.0./1.0e+2,0.0,0.0,1.0,0.0,t18.*t22.*(-8.333333333333333e-8)-(t22.*t33)./1.2e+7,t19.*t22.*(-8.333333333333333e-8)-(t22.*t35)./1.2e+7];
mt2 = [1.0./1.0e+2,t38+t51-(t22.*t32)./6.0e+4,t39+t52-(t22.*t34)./6.0e+4,1.0];
fx = reshape([mt1,mt2],6,6);
if nargout > 1
    fu = reshape([t23+t29+t48+t57,t20+t27+t46+t58,-5.208333333333333e-4,t14+t26+t44+t45+t53+t54,t8+t24+t40+t41+t55+t56,-5.0./4.8e+1,t23+t29+t48-t57,t20+t27+t46-t58,5.208333333333333e-4,t14+t26+t44+t45-t53-t54,t8+t24+t40+t41-t55-t56,5.0./4.8e+1],[6,2]);
end
