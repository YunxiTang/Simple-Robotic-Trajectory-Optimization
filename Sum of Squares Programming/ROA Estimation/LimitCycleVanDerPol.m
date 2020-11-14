%%% checking the ROA by repeatedly forward simulations
clear all;
tspan = [0,60];
n = 150;
q1 = linspace(-1.5,1.5,n);
N = n^2;
q2 = q1;
count = 0;
figure(1);

for i=1:length(q1)
    for j=1:length(q2)
        count = count + 1;
        x_ini = [q1(i);q2(j)];
        SOL = ode45(@(t,x)(VanderDyna(t,x)),tspan,x_ini);
        x = SOL.y;
        if sqrt(x(1,end)^2+x(2,end)^2) <= 1e-3
            plot([q1(i)],[q2(j)],'Color','yellow','Marker','.','LineWidth',0.1);
            disp('Point In.');
%             plot(x(1,:),x(2,:),'r');
        else
            h=plot([q1(i)],[q2(j)],'Color','green','Marker','.','LineWidth',0.1);            
            disp('Point Out.');
        end
       
        hold on;
        fprintf('Progress is %d/%d now. \n', count,N);
    end
end

xlim([-1.6,1.6]);
ylim([-1.6,1.6]);
axis equal;
grid on;
xlabel('x1');ylabel('x2');
title('ROA');
alpha(0.01);
saveas(gcf, 'ROA', 'fig')
%%
