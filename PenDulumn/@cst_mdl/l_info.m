function [l,lx,lu,lxx,lux,lxu,luu] = l_info(in1,u1)
%L_INFO
%    [L,LX,LU,LXX,LUX,LXU,LUU] = L_INFO(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    22-Mar-2021 09:45:33

x1 = in1(1,:);
x2 = in1(2,:);
l = ((x1./2.0-1.57e+2./1.0e+2).*(x1-1.57e+2./5.0e+1))./1.0e+2+u1.^2./2.0e+2+x2.^2./2.0e+2;
if nargout > 1
    lx = [x1./1.0e+2-3.14e-2;x2./1.0e+2];
end
if nargout > 2
    lu = u1./1.0e+2;
end
if nargout > 3
    lxx = reshape([1.0./1.0e+2,0.0,0.0,1.0./1.0e+2],[2,2]);
end
if nargout > 4
    lux = [0.0,0.0];
end
if nargout > 5
    lxu = [0.0;0.0];
end
if nargout > 6
    luu = 1.0./1.0e+2;
end
