function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    13-May-2021 18:31:45

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = x2-3.0;
t3 = t2.*1.25e+2;
lf = t2.*t3+(x1.*1.25e+2-6.25e+2./2.0).*(x1-5.0./2.0)+(x3-pi./2.0).*(x3.*1.25e+2-pi.*(1.25e+2./2.0))+x4.^2.*1.25e+2;
if nargout > 1
    lfx = [x1.*2.5e+2-6.25e+2;t2.*2.5e+2;x3.*2.5e+2-pi.*1.25e+2;x4.*2.5e+2];
end
if nargout > 2
    lfxx = reshape([2.5e+2,0.0,0.0,0.0,0.0,2.5e+2,0.0,0.0,0.0,0.0,2.5e+2,0.0,0.0,0.0,0.0,2.5e+2],[4,4]);
end
