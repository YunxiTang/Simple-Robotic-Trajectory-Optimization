function final_c = final_c(in1)
%FINAL_C
%    FINAL_C = FINAL_C(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    08-Jun-2021 16:02:23

x1 = in1(1,:);
x2 = in1(2,:);
t2 = pi./4.5e+1;
t3 = -t2;
final_c = [t3+x1;t2-x1;t2+x2;t3-x2];
