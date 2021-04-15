function [fx,fu] = getLinSys(in1,in2)
%GETLINSYS
%    [FX,FU] = GETLINSYS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    13-Apr-2021 19:02:45

u1 = in2(1,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = cos(x3);
t3 = sin(x3);
fx = reshape([1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,(t2.*x4)./1.0e+2,t3.*x4.*(-1.0./1.0e+2),1.0,0.0,t3./1.0e+2,t2./1.0e+2,u1./1.0e+2,1.0],[4,4]);
if nargout > 1
    fu = reshape([0.0,0.0,x4./1.0e+2,0.0,0.0,0.0,0.0,1.0./1.0e+2],[4,2]);
end
