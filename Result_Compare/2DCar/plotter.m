clear;
close all;
clc;
%% time
% time = load('time.mat');
% t = time.t;
% t_ss = linspace(0, 3, 150);
%% al ms ddp
x_almsddp = load('x_al.mat');
u_almsddp = load('u_al.mat');

x_al = x_almsddp.xsol;
u_al = u_almsddp.usol;
%% rlb ms ddp
x_rlbmsddp = load('x_rlb.mat');
u_rlbmsddp = load('u_rlb.mat');
x_rlb = x_rlbmsddp.xsol;
u_rlb = u_rlbmsddp.usol;

%% single sshooting
x_sshooting = load('x_ss.mat');
u_sshooting = load('u_ss.mat');

x_ss = x_sshooting.xsol;
u_ss = u_sshooting.u_opt;

%% dircol 
x_dircol = load('x_dircol.mat');
u_dircol = load('u_dircol.mat');

x_dir = x_dircol.z;
u_dir = u_dircol.u;

%% single rlb
x_Srlb = load('x_srlb.mat');
u_Srlb = load('u_srlb.mat');
x_srlb = x_Srlb.xsol;
u_srlb = u_Srlb.usol;

%%
Obstacles = [0.0 1.3 2.0;
             1.0 1.3 2.5;
             0.5 0.6 0.4];
%% 
figure(1);
plot_obstacle(Obstacles, 1);
h1=plot(x_al(1,1), x_al(2,1), 'gd', 'LineWidth', 5.0); hold on;
h2=plot(2.5, 3.0, 'rx', 'LineWidth', 12.0); hold on;
h3=plot(x_al(1,:), x_al(2,:), 'k-.','LineWidth',3.5); hold on;
h4=plot(x_srlb(1,:), x_srlb(2,:), 'Color','#4DBEEE', 'LineStyle', '-.','LineWidth',3.0);
h5=plot(x_ss(1,:), x_ss(2,:), 'b-.', 'LineWidth',3.5); hold on;
h6=plot(x_dir(1,:), x_dir(2,:), 'Color','#7E2F8E', 'LineStyle', '-.','LineWidth',3.5);hold on;
h7=plot(x_rlb(1,:), x_rlb(2,:), 'r-','LineWidth',3.0); hold on;
axis equal;
xlabel('$x[m]$', 'Interpreter', 'latex', 'FontSize', 18, 'FontAngle','italic');
ylabel('$y[m]$', 'Interpreter', 'latex', 'FontSize', 18, 'FontAngle','italic');
title('$2D\;Car$', 'Interpreter', 'latex', 'FontSize', 20, 'FontAngle','italic');
ld2 = legend([h1,h2,h3,h4,h5,h6,h7],'$\textbf{x}_0$','$\textbf{x}_f$','$AL$', '$SRLB$','$DSS$','$DIRCOL$', '$Ours$',...
    'Interpreter','latex', 'FontSize', 12,'NumColumns',2);
set (gcf,'Position',[400,100,500,500], 'color','w');
% 
% figure(2);
% plot(t(1:end-1), u_al, 'k-.','LineWidth',2.0);hold on;
% plot(t(1:end-1), u_rlb, 'r-','LineWidth',2.0); hold on;
% plot(t_ss, u_ss, 'b-.','LineWidth',2.0); hold on;
% plot(t(1:end-1), u_dir, 'Color','#7E2F8E', 'LineStyle', '-.','LineWidth',2.0); hold on;
% plot(t(1:end-1), u_srlb, 'Color','#4DBEEE', 'LineStyle', '-.','LineWidth',2.0); 
% xlabel('$t[s]$', 'Interpreter', 'latex', 'FontSize', 18, 'FontAngle','italic');
% ylabel('$u_F[N]$', 'Interpreter', 'latex', 'FontSize', 18, 'FontAngle','italic');
% title('$CartPole$', 'Interpreter', 'latex', 'FontSize', 20, 'FontAngle','italic');
% legend('$AL-DDP$', '$Ours$','$SS$','$DIRCOL$', '$SRLB$', 'Interpreter','latex', 'FontSize', 15);
% set (gcf,'Position',[400,100,600,500], 'color','w');