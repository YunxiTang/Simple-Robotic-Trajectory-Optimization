function final_c = final_c(in1)
%FINAL_C
%    FINAL_C = FINAL_C(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    22-Apr-2021 14:39:11

x1 = in1(1,:);
x2 = in1(2,:);
t2 = x2-pi;
t3 = -t2;
final_c = [x1;-x1;t2;t3];
