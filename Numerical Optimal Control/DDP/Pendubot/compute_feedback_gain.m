Quu = sol.Q_uu;
Qux = sol.Q_ux;
Qu = sol.Q_u;
N = length(sol.t);
M = N-1;
K = zeros(M,4);

Q_uu = zeros(M,1);
Q_u = zeros(M,1);
Q_ux = zeros(M,4);
for i = 1:M
    Q_uu(i,:) = Quu{i}(1);
    Q_u(i,:) = Qu{i}(1);
    Q_ux(i,:) = Qux{i};
    k = -inv(Quu{i}) * Qux{i};
    K(i,:) = -k;
end
figure;
plot(sol.t(1:M),K,'LineWidth',2.0);
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
grid on;
title('Feedback Gains','FontName',"latex", "FontSize", 20);
figure;
plot(sol.t(1:M),Q_uu,'LineWidth',2.0);
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
grid on;
title('Q_{uu}','FontName',"latex", "FontSize", 20);

figure;
plot(sol.t(1:M),Q_u,'LineWidth',2.0);
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
grid on;
title('Q_{u}','FontName',"latex", "FontSize", 20);

figure;
plot(sol.t(1:M),Q_ux,'LineWidth',2.0);
xlabel("Time [s]", "Interpreter", "latex", "FontSize", 20);
grid on;
title('Q_{ux}','FontName',"latex", "FontSize", 20);