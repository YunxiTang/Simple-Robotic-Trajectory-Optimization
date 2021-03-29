function [l,lx,lu,lxx,lux,lxu,luu] = l_info(in1,in2)
%L_INFO
%    [L,LX,LU,LXX,LUX,LXU,LUU] = L_INFO(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    26-Mar-2021 14:49:26

u1 = in2(1,:);
u2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
l = u1.^2./1.0e+3+u2.^2./1.0e+3+x1.^2./1.0e+3+x2.^2./1.0e+3+x3.^2./1.0e+3+x4.^2./1.0e+3+x5.^2./1.0e+3+x6.^2./1.0e+3;
if nargout > 1
    lx = [x1./5.0e+2;x2./5.0e+2;x3./5.0e+2;x4./5.0e+2;x5./5.0e+2;x6./5.0e+2];
end
if nargout > 2
    lu = [u1./5.0e+2;u2./5.0e+2];
end
if nargout > 3
    lxx = reshape([1.0./5.0e+2,0.0,0.0,0.0,0.0,0.0,0.0,1.0./5.0e+2,0.0,0.0,0.0,0.0,0.0,0.0,1.0./5.0e+2,0.0,0.0,0.0,0.0,0.0,0.0,1.0./5.0e+2,0.0,0.0,0.0,0.0,0.0,0.0,1.0./5.0e+2,0.0,0.0,0.0,0.0,0.0,0.0,1.0./5.0e+2],[6,6]);
end
if nargout > 4
    lux = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[2,6]);
end
if nargout > 5
    lxu = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[6,2]);
end
if nargout > 6
    luu = reshape([1.0./5.0e+2,0.0,0.0,1.0./5.0e+2],[2,2]);
end
