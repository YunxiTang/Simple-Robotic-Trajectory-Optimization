function [cx,cu] = algrad(in1,in2)
%ALGRAD
%    [CX,CU] = ALGRAD(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    18-May-2021 15:30:23

x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1.*2.0;
t3 = x2.*2.0;
t4 = -t2;
t5 = -t3;
cx = reshape([t4,t4,t4,t5,t5,t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,6]);
if nargout > 1
    cu = reshape([0.0,0.0,0.0,0.0,0.0,0.0],[3,2]);
end
