%%% Test...ing
clear all;
clc;

%% YAMLIP
sdpvar x y lower
p = (1+x*y)^2-x*y+(1-y)^2;
F = sos(p-lower);
solvesos(F,-lower,[],lower);
value(lower)