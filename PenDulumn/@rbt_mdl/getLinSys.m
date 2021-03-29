function [fx,fu,fx_ctn,fu_ctn] = getLinSys(in1,u1)
%GETLINSYS
%    [FX,FU,FX_CTN,FU_CTN] = GETLINSYS(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    26-Mar-2021 11:16:52

x1 = in1(1,:);
t2 = cos(x1);
fx = reshape([1.0,t2.*(-8.456896551724138e-2),1.0./1.0e+2,9.979310344827586e-1],[2,2]);
if nargout > 1
    fu = [0.0;1.0./5.8e+1];
end
if nargout > 2
    fx_ctn = reshape([0.0,t2.*(-9.81e+2./1.16e+2),1.0,-6.0./2.9e+1],[2,2]);
end
if nargout > 3
    fu_ctn = [0.0;5.0e+1./2.9e+1];
end
