function final_c = final_c(in1)
%FINAL_C
%    FINAL_C = FINAL_C(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    17-May-2021 18:50:36

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t2 = x1-1.0;
t3 = x2-3.0./2.0;
final_c = [t2;-t2;t3;-t3;x3;-x3];
