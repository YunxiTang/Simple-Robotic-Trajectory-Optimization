%% define a structure
clc;clear all;
sdpvar x y

problem.Lyafunc = x^2 + y^2;
problem.Controller = polynomial(x,1);