function c = c(in1,in2)
%C
%    C = C(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    13-May-2021 11:37:18

u1 = in2(1,:);
u2 = in2(2,:);
x1 = in1(1,:);
x2 = in1(2,:);
c = [-(x2-1.0).^2-x1.^2+1.0./4.0;-(x1-1.3e+1./1.0e+1).^2-(x2-1.3e+1./1.0e+1).^2+9.0./2.5e+1;-(x1-2.0).^2-(x2-5.0./2.0).^2+1.0./4.0;u1-9.0./2.0;-u1-9.0./2.0;u2-9.0./2.0;-u2-9.0./2.0];
