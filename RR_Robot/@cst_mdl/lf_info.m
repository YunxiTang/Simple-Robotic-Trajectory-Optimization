function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    08-Apr-2021 21:34:40

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = x1-pi;
t3 = t2.*6.0;
lf = t2.*t3+x2.^2.*6.0+x3.^2.*6.0+x4.^2.*6.0;
if nargout > 1
    lfx = [t2.*1.2e+1,x2.*1.2e+1,x3.*1.2e+1,x4.*1.2e+1];
end
if nargout > 2
    lfxx = reshape([1.2e+1,0.0,0.0,0.0,0.0,1.2e+1,0.0,0.0,0.0,0.0,1.2e+1,0.0,0.0,0.0,0.0,1.2e+1],[4,4]);
end
