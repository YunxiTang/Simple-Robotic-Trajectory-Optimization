function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    23-Mar-2021 18:33:04

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
t2 = x2-1.0;
t3 = t2.*2.5e+1;
lf = t2.*t3+(x1.*2.5e+1-4.75e+2./2.0).*(x1-1.9e+1./2.0)+x3.^2.*2.5e+1+x4.^2.*2.5e+1+x5.^2.*2.5e+1+x6.^2.*2.5e+1;
if nargout > 1
    lfx = [x1.*5.0e+1-4.75e+2,t2.*5.0e+1,x3.*5.0e+1,x4.*5.0e+1,x5.*5.0e+1,x6.*5.0e+1];
end
if nargout > 2
    lfxx = reshape([5.0e+1,0.0,0.0,0.0,0.0,0.0,0.0,5.0e+1,0.0,0.0,0.0,0.0,0.0,0.0,5.0e+1,0.0,0.0,0.0,0.0,0.0,0.0,5.0e+1,0.0,0.0,0.0,0.0,0.0,0.0,5.0e+1,0.0,0.0,0.0,0.0,0.0,0.0,5.0e+1],[6,6]);
end
