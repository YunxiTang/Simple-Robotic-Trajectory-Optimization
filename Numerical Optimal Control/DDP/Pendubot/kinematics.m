function [p1,p2] = kinematics(l1, l2, x)
%KINEMATICS forward kinematics

q1 = x(1,:);
q2 = x(2,:);

p1 = [l1*sin(q1);
      -l1*cos(q1)];
p2 = p1 + [l2*sin(q1+q2);
           -l2*cos(q1+q2)];

end
