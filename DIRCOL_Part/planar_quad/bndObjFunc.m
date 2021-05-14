function [J] = bndObjFunc(t0,x0,tf,xf, params)
%BNDOBJFUNC 
% boundary objective function
Jf = (xf - params.xf)' * params.Qf * (xf - params.xf) / 2;
J0 = 0;%(x0 - params.xf)' * params.Q * (x0 - params.xf) / 2;
J = J0+ Jf;
end

