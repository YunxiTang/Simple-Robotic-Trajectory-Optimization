function cx = cx(in1,in2)
%CX
%    CX = CX(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    13-May-2021 11:37:18

x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1.*2.0;
t3 = x2.*2.0;
t4 = -t2;
t5 = -t3;
cx = reshape([t4,t4+1.3e+1./5.0,t4+4.0,0.0,0.0,0.0,0.0,t5+2.0,t5+1.3e+1./5.0,t5+5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[7,4]);
