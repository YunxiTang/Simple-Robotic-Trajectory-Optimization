function [lf,lfx,lfxx] = lf_info(in1)
%LF_INFO
%    [LF,LFX,LFXX] = LF_INFO(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    22-Mar-2021 09:45:33

x1 = in1(1,:);
x2 = in1(2,:);
lf = (x1.*5.0-1.57e+2./1.0e+1).*(x1-1.57e+2./5.0e+1)+x2.^2.*5.0;
if nargout > 1
    lfx = [x1.*1.0e+1-1.57e+2./5.0,x2.*1.0e+1];
end
if nargout > 2
    lfxx = reshape([1.0e+1,0.0,0.0,1.0e+1],[2,2]);
end
