function [lf,lfx,lfxx] = lf_info(in1,in2)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    18-May-2021 21:15:13

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
xref1 = in2(1,:);
xref2 = in2(2,:);
xref3 = in2(3,:);
xref4 = in2(4,:);
xref5 = in2(5,:);
xref6 = in2(6,:);
t2 = -xref6;
t3 = t2+x6;
lf = (x3.*5.0-xref3.*5.0).*(x3-xref3)+(x4.*2.5e+1-xref4.*2.5e+1).*(x4-xref4)+(x5.*2.5e+1-xref5.*2.5e+1).*(x5-xref5)+(x1.*6.0e+1-xref1.*6.0e+1).*(x1-xref1)+(x2.*6.0e+1-xref2.*6.0e+1).*(x2-xref2)+t3.*(x6./2.0-xref6./2.0);
if nargout > 1
    lfx = [x1.*1.2e+2-xref1.*1.2e+2,x2.*1.2e+2-xref2.*1.2e+2,x3.*1.0e+1-xref3.*1.0e+1,x4.*5.0e+1-xref4.*5.0e+1,x5.*5.0e+1-xref5.*5.0e+1,t3];
end
if nargout > 2
    lfxx = reshape([1.2e+2,0.0,0.0,0.0,0.0,0.0,0.0,1.2e+2,0.0,0.0,0.0,0.0,0.0,0.0,1.0e+1,0.0,0.0,0.0,0.0,0.0,0.0,5.0e+1,0.0,0.0,0.0,0.0,0.0,0.0,5.0e+1,0.0,0.0,0.0,0.0,0.0,0.0,1.0],[6,6]);
end
