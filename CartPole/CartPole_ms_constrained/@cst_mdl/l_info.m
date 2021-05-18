function [l,lx,lu,lxx,lux,lxu,luu] = l_info(in1,u1)
%L_INFO
%    [L,LX,LU,LXX,LUX,LXU,LUU] = L_INFO(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    17-May-2021 10:24:09

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
l = ((x2-pi).*(x2./2.0e+4-pi./2.0e+4))./1.0e+2+u1.^2./2.0e+3+x1.^2./2.0e+6+x3.^2./2.0e+6+x4.^2./2.0e+6;
if nargout > 1
    lx = [x1./1.0e+6;x2./1.0e+6-pi./1.0e+6;x3./1.0e+6;x4./1.0e+6];
end
if nargout > 2
    lu = u1./1.0e+3;
end
if nargout > 3
    lxx = reshape([1.0e-6,0.0,0.0,0.0,0.0,1.0e-6,0.0,0.0,0.0,0.0,1.0e-6,0.0,0.0,0.0,0.0,1.0e-6],[4,4]);
end
if nargout > 4
    lux = [0.0,0.0,0.0,0.0];
end
if nargout > 5
    lxu = [0.0;0.0;0.0;0.0];
end
if nargout > 6
    luu = 1.0./1.0e+3;
end
