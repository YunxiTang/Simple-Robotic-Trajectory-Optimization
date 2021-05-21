function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    21-May-2021 16:18:53

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = x2-pi;
lf = t2.*(x2.*(5.0./2.0)-pi.*(5.0./2.0))+x1.^2.*(5.0./2.0)+x3.^2.*(5.0./2.0)+x4.^2.*(5.0./2.0);
if nargout > 1
    lfx = [x1.*5.0;t2.*5.0;x3.*5.0;x4.*5.0];
end
if nargout > 2
    lfxx = reshape([5.0,0.0,0.0,0.0,0.0,5.0,0.0,0.0,0.0,0.0,5.0,0.0,0.0,0.0,0.0,5.0],[4,4]);
end
