clear;
clc;
disp('~~~~~~~~~~!Testing EoM!~~~~~~~~~');

q_test = [0.9134    0.6324    0.0975];
dq_test = [0.2785    0.5469    0.9575];
x_test = [q_test dq_test]';
C_test = [0  -0.132706376041645   1.037364687878396
   0.067578580595352                   0                   0
  -0.301729572401184                   0                   0];

M_test = [ 6.437500000000000  -0.840681276908014   1.019254579560482
          -0.840681276908014   0.437500000000000                   0
          1.019254579560482                   0   0.520625000000000];
G_test = 100 * [-1.067750441505868
            0.101474058206939
            -0.028410069075374];
   
[M, C, G, B] = EoM(x_test);

error_C = round(C - C_test, 15);
fprintf('error_C: \n');
disp(error_C)


error_M = round(M - M_test, 15);
fprintf('error_M: \n');
disp(error_M)

error_G = round(G - G_test, 15);
fprintf('error_G: \n');
disp(error_G)

disp('~~~~~~~~~~!Testing Done!~~~~~~~~~');