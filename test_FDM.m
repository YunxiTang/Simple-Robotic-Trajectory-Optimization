%%%%% Finite Difference Method to linearize system dynamics
x0 = [1;1];
J0 =  Finite_Diff(@(x)F(x),x0);
N = 500;
X = linspace(0.0, pi, N);
Y = linspace(0.0, pi, N);
Q = [X;
     Y];
FDM_J = zeros(4, N);
Ana_J = zeros(4, N);
for i=1:N
    q = [X(i);Y(i)];
    J_fdm = Finite_Diff(@(x)F(x),q,1e-3);
    J_ana = Analytical_Jacob(q);
    FDM_J(:,i) = reshape(J_fdm,4,1);
    Ana_J(:,i) = reshape(J_ana,4,1);
end
figure();
yyaxis left;h2=plot(1:N,Ana_J,'r-','LineWidth',4.0);hold on;
h1=plot(1:N,FDM_J,'b-','LineWidth',1.5);

grid on;

xlabel('$Sample\;Points$',...
       'Interpreter','latex','Fontsize',14);
ylabel('$Jacobians$','Interpreter','latex','Fontsize',14);

yyaxis right;
h3=plot(1:N,FDM_J-Ana_J,'k--','LineWidth',1.0);
legend([h1(1),h2(1),h3(1)],'$FDM\;Jacobian$','$Ground\;Truth$','$Error$',...
       'Interpreter','latex','Fontsize',14);
%% functions
function [dq] = F(q)
    x1 = q(1);
    x2 = q(2);
    dq1 = sin(x1)   + cos(x2);
    dq2 = cos(2*x1) + sin(3*x2);
    dq = [dq1;dq2];
end

function J = Finite_Diff(func, qbar, dq)
    if nargin == 3
        h = dq;
    else
        h = 1e-4;
    end
    Nx = numel(qbar);
    J = zeros(Nx,Nx);
    H = eye(Nx) * h;
    for i=1:Nx
        J(:,i) = (func(qbar + H(:,i)) - func(qbar - H(:,i))) / (2*h);
    end
end

function J = Analytical_Jacob(q)
    x1 = q(1);
    x2 = q(2);
    J = [cos(x1)        -sin(x2);
         -2*sin(2*x1) 3*cos(3*x2)];
end