function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    07-Apr-2021 19:48:23

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
lf = (x1.*6.0e+1-9.42e+2./5.0).*(x1-1.57e+2./5.0e+1)+x2.^2.*6.0e+1+x3.^2.*6.0e+1+x4.^2.*6.0e+1;
if nargout > 1
    lfx = [x1.*1.2e+2-3.768e+2,x2.*1.2e+2,x3.*1.2e+2,x4.*1.2e+2];
end
if nargout > 2
    lfxx = reshape([1.2e+2,0.0,0.0,0.0,0.0,1.2e+2,0.0,0.0,0.0,0.0,1.2e+2,0.0,0.0,0.0,0.0,1.2e+2],[4,4]);
end
