function final_c = final_c(in1)
%FINAL_C
%    FINAL_C = FINAL_C(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    16-Apr-2021 15:41:59

x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1-3.0;
t3 = x2-3.0;
final_c = [t2;-t2;t3;-t3];
