function c = c(in1,u1)
%C
%    C = C(IN1,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    06-May-2021 11:18:23

x1 = in1(1,:);
x2 = in1(2,:);
t2 = pi.*(3.0./2.0);
t3 = -t2;
c = [u1-2.5e+1;-u1-2.5e+1;t3+x2;t3-x2;x1-4.0./5.0;-x1-4.0./5.0];
