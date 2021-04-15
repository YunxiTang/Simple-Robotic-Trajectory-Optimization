function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    15-Apr-2021 15:24:04

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x6 = in1(6,:);
x7 = in1(7,:);
x8 = in1(8,:);
lf = x1.^2.*3.0+x2.^2.*3.0+x3.^2.*3.0+x6.^2.*3.0+x7.^2.*3.0+x8.^2.*3.0;
if nargout > 1
    lfx = [x1.*6.0;x2.*6.0;x3.*6.0;0.0;0.0;x6.*6.0;x7.*6.0;x8.*6.0;0.0;0.0];
end
if nargout > 2
    lfxx = reshape([6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,6.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[10,10]);
end
