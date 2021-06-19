clc;
clear;
close all;
%% Hyb-CDDP: x=0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7 (M=5)
hyb_cddp_s1_iter = [10     5      5      5      5      5      5];
hyb_cddp_s1_time = [1.6924 1.8811 1.9737 1.7688 1.7777 1.7659 1.7720];
hyb_cddp_s1_ave  = hyb_cddp_s1_time ./ hyb_cddp_s1_iter;

hyb_cddp_s2_iter = [149    122    53     22     51     117    24];
hyb_cddp_s2_time = [6.2326 6.5245 2.1220 1.0214 2.1785 4.8376 1.1482];
hyb_cddp_s2_ave  = hyb_cddp_s2_time ./ hyb_cddp_s2_iter;

hyb_cddp_iter    = hyb_cddp_s1_iter + hyb_cddp_s2_iter;
hyb_cddp_time    = hyb_cddp_s1_time + hyb_cddp_s2_time;
hyb_cddp_ave     = hyb_cddp_time ./ hyb_cddp_iter;

plot(hyb_cddp_s1_ave);hold on;
plot(hyb_cddp_s2_ave);hold on;
plot(hyb_cddp_ave);

figure(2);
bar(hyb_cddp_s1_ave);hold on;
bar(hyb_cddp_ave);hold on;
bar(-hyb_cddp_s2_ave);
