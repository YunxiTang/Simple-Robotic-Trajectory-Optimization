function [fx,fu,fx_ctn,fu_ctn] = getLinSys(in1,u1)
%GETLINSYS
%    [FX,FU,FX_CTN,FU_CTN] = GETLINSYS(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    22-Apr-2021 17:29:00

x1 = in1(1,:);
t2 = cos(x1);
fx = reshape([1.0,t2.*(-5.637931034482759e-2),1.0./1.0e+2,7.24e+2./7.25e+2],[2,2]);
if nargout > 1
    fu = [0.0;1.0./8.7e+1];
end
if nargout > 2
    fx_ctn = reshape([0.0,t2.*(-3.27e+2./5.8e+1),1.0,-4.0./2.9e+1],[2,2]);
end
if nargout > 3
    fu_ctn = [0.0;1.0e+2./8.7e+1];
end
