function [M, C, G, B] = EoM(x)
%EOM 
    q1 = x(1);
    q2 = x(2);
    q3 = x(3);
    dq1 = x(4);
    dq2 = x(5);
    dq3 = x(6);
    
    [m1, m2, m3, l1, l2, l3, g] = set_parameters();
    M = eval_M_tmp(l1,l2,l3,m1,m2,m3,q1,q2,q3);
    C = eval_C_tmp(dq1,dq2,dq3,l1,l2,l3,m2,m3,q1,q2,q3);
    G = eval_G_tmp(g,l1,l2,l3,m1,m2,m3,q1,q2,q3);
    B = eval_B_tmp();
end

