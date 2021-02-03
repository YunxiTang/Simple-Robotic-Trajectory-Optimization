function [] = plot_trail(time_seq,trail,defect,til)
%PLOT_TRAIL plot one trial WITH defect
    figure;
    defect  =[[0 0 0 0]; defect];
    [S, ~, Nx] = size(trail);
    for k=1:S 
        for i=1:Nx
            subplot(Nx,1,i);
            plot(squeeze(time_seq(k,:,1)),squeeze(trail(k,:,i)),'k-','LineWidth',2.0);hold on;
            plot(squeeze(time_seq(k,1,1)),defect(k,i,1),'o',...
                'MarkerEdgeColor','b','MarkerFaceColor','r','LineWidth',1.0);
            hold on;
            grid on;
            ylabel(strcat('$x_',num2str(i),'$'),'Interpreter','latex','Fontsize',12);
        end
    end
    xlabel('$Time \; [s]$','Interpreter','latex','Fontsize',12);
    subplot(Nx,1,1);
    title(til,'Interpreter','latex','Fontsize',15);
    hold off;
    drawnow;
end

