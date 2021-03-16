function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    16-Mar-2021 16:51:49

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = x1.*2.5e+3-6.25e+3;
t3 = x2.*2.5e+3-6.25e+3;
lf = t2.*(x1-5.0./2.0)+t3.*(x2-5.0./2.0)+x3.^2.*2.5e+3+x4.^2.*2.5e+3;
if nargout > 1
    lfx = [t2.*2.0;t3.*2.0;x3.*5.0e+3;x4.*5.0e+3];
end
if nargout > 2
    lfxx = reshape([5.0e+3,0.0,0.0,0.0,0.0,5.0e+3,0.0,0.0,0.0,0.0,5.0e+3,0.0,0.0,0.0,0.0,5.0e+3],[4,4]);
end
