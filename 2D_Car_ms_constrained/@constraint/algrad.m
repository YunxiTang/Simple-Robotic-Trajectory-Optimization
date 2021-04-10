function [cx,cu] = algrad(in1,in2)
%ALGRAD
%    [CX,CU] = ALGRAD(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    10-Apr-2021 17:58:55

x1 = in1(1,:);
x2 = in1(2,:);
cx = [x1.*-2.0+3.0,x2.*-2.0+3.0,0.0,0.0];
if nargout > 1
    cu = [0.0,0.0];
end
