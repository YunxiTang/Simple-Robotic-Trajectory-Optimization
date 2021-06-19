function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    28-May-2021 14:00:23

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
t2 = x1-2.0;
t3 = x2-1.0;
lf = t2.*(x1.*(5.0./2.0)-5.0)+t3.*(x2.*(5.0./2.0)-5.0./2.0)+x3.^2.*(5.0./2.0)+x4.^2.*(5.0./2.0)+x5.^2.*(5.0./2.0)+x6.^2.*(5.0./2.0);
if nargout > 1
    lfx = [t2.*5.0;t3.*5.0;x3.*5.0;x4.*5.0;x5.*5.0;x6.*5.0];
end
if nargout > 2
    lfxx = reshape([5.0,0.0,0.0,0.0,0.0,0.0,0.0,5.0,0.0,0.0,0.0,0.0,0.0,0.0,5.0,0.0,0.0,0.0,0.0,0.0,0.0,5.0,0.0,0.0,0.0,0.0,0.0,0.0,5.0,0.0,0.0,0.0,0.0,0.0,0.0,5.0],[6,6]);
end
