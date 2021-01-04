function runAnimation(t, p1, p2)

    %ANIMATION 
    figure(3333);
    for i=1:2:length(t)
        clf;
        plot([0 p1(1,i)],[0 p1(2,i)],'Color','k','Linewidth',4); % link1
        hold on;

        plot([p1(1,i) p2(1,i)],[p1(2,i) p2(2,i)],'Color','r','Linewidth',3); % link2
        hold on;
        scatter(p1(1,i), p1(2,i),60,'b', 'filled');
        hold on;
        scatter([0],[0],'Linewidth',3);
        hold on;
        axis equal;

        xlim([-2,2]);
        ylim([-2,2]);
        pause(0.025);
        grid on;
    end
    plot(p1(1,:),p1(2,:),'.-','linewidth',0.5); % trajectory 1
    hold on;
    plot(p2(1,:),p2(2,:),'.-'); % trajectory 2
    axis equal;
    title('Pendubot','Interpreter','latex','FontSize',15);
end