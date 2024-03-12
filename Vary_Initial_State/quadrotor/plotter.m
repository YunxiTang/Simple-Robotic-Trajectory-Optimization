clear;
close all;
clc;

%% al ms ddp
x_almsddp = load('x_al_1.mat');
x_al_1 = x_almsddp.xsol;

x_almsddp = load('x_al_2.mat');
x_al_2 = x_almsddp.xsol;

x_almsddp = load('x_al_3.mat');
x_al_3 = x_almsddp.xsol;

x_almsddp = load('x_al_4.mat');
x_al_4 = x_almsddp.xsol;

x_almsddp = load('x_al_5.mat');
x_al_5 = x_almsddp.xsol;

x_almsddp = load('x_al_6.mat');
x_al_6 = x_almsddp.xsol;

x_almsddp = load('x_al_7.mat');
x_al_7 = x_almsddp.xsol;
%% rlb ms ddp
x_rlbmsddp = load('x_rlb_1.mat');
x_rlb_1 = x_rlbmsddp.xsol;

x_rlbmsddp = load('x_rlb_2.mat');
x_rlb_2 = x_rlbmsddp.xsol;

x_rlbmsddp = load('x_rlb_3.mat');
x_rlb_3 = x_rlbmsddp.xsol;

x_rlbmsddp = load('x_rlb_4.mat');
x_rlb_4 = x_rlbmsddp.xsol;

x_rlbmsddp = load('x_rlb_5.mat');
x_rlb_5 = x_rlbmsddp.xsol;

x_rlbmsddp = load('x_rlb_6.mat');
x_rlb_6 = x_rlbmsddp.xsol;

x_rlbmsddp = load('x_rlb_7.mat');
x_rlb_7 = x_rlbmsddp.xsol;
%%
Obstacles = [2.0 3.0 4.0;
             2.0 1.0 2.0;
             0.6 0.5 0.6];
%% 
figure(1);
plot_obstacle(Obstacles, 1);
%%
h1=plot(1.0, 1.5, 'gx', 'LineWidth', 5.0); hold on;
%%
h2=plot(x_al_1(1,1), x_al_1(2,1), 'mx', 'LineWidth', 5.0); hold on;
h3=plot(x_al_2(1,1), x_al_2(2,1), 'mx', 'LineWidth', 5.0); hold on;
h4=plot(x_al_3(1,1), x_al_3(2,1), 'mx', 'LineWidth', 5.0); hold on;
h5=plot(x_al_4(1,1), x_al_4(2,1), 'mx', 'LineWidth', 5.0); hold on;
h6=plot(x_al_5(1,1), x_al_5(2,1), 'mx', 'LineWidth', 5.0); hold on;
h7=plot(x_al_6(1,1), x_al_6(2,1), 'mx', 'LineWidth', 5.0); hold on;
h8=plot(x_al_7(1,1), x_al_7(2,1), 'mx', 'LineWidth', 5.0); hold on;
%% case 1
pk1=plot(x_al_1(1,:), x_al_1(2,:), 'r-.','LineWidth',2.0); hold on;
plot(x_rlb_1(1,:), x_rlb_1(2,:), 'r-','LineWidth',2.0); hold on;
%% case 2
pk2=plot(x_al_2(1,:), x_al_2(2,:), 'b-.','LineWidth',2.0); hold on;
plot(x_rlb_2(1,:), x_rlb_2(2,:), 'b-','LineWidth',2.0); hold on;
%% case 3
pk3=plot(x_al_3(1,:), x_al_3(2,:), ' k-.','LineWidth',2.0); hold on;
plot(x_rlb_3(1,:), x_rlb_3(2,:), 'k-','LineWidth',2.0); hold on;
%% case 3
pk4=plot(x_al_3(1,:), x_al_3(2,:), ' k-.','LineWidth',2.0); hold on;
plot(x_rlb_3(1,:), x_rlb_3(2,:), 'k-','LineWidth',2.0); hold on;
%% case 4
pk5=plot(x_al_4(1,:), x_al_4(2,:), 'Color', '#0072BD', 'LineStyle','-.','LineWidth',2.0); hold on;
plot(x_rlb_4(1,:), x_rlb_4(2,:), 'Color', '#0072BD','LineWidth',2.0); hold on;

%% case 5
pk6=plot(x_al_5(1,:), x_al_5(2,:), 'c-.','LineWidth',2.0); hold on;
plot(x_rlb_5(1,:), x_rlb_5(2,:), 'c-','LineWidth',2.0); hold on;

%% case 6
pk7=plot(x_al_6(1,:), x_al_6(2,:), 'Color', '#7E2F8E', 'LineStyle','-.','LineWidth',2.0); hold on;
plot(x_rlb_6(1,:), x_rlb_6(2,:), 'Color', '#7E2F8E','LineWidth',2.0); hold on;

%% case 6
pk8=plot(x_al_7(1,:), x_al_7(2,:), 'Color', '#77AC30', 'LineStyle','-.','LineWidth',2.0); hold on;
plot(x_rlb_7(1,:), x_rlb_7(2,:), 'Color', '#77AC30','LineWidth',2.0); hold on;
%%
vis = 0.5;
pk1.Color(4)=vis;pk2.Color(4)=vis;
pk3.Color(4)=vis;pk4.Color(4)=vis;
pk5.Color(4)=vis;pk6.Color(4)=vis;
pk7.Color(4)=vis;pk8.Color(4)=vis;
%%
axis equal;
xlabel('$x[m]$', 'Interpreter', 'latex', 'FontSize', 18, 'FontAngle','italic');
ylabel('$y[m]$', 'Interpreter', 'latex', 'FontSize', 18, 'FontAngle','italic');
title('$Planar\;Quadrotor$', 'Interpreter', 'latex', 'FontSize', 20, 'FontAngle','italic');
% ld2 = legend([h1,h2,h3,h4],'$\textbf{x}_f$','$\textbf{x}_{01}$','$\textbf{x}_{02}$','$\textbf{x}_{03}$',...
%     'Interpreter','latex', 'FontSize', 15,'NumColumns',2);
set (gcf,'Position',[400,100,500,500], 'color','w');
