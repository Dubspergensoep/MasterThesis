function parforplot_a2f2color(rho,timecell,fcell,dS,sig,plottitle,saveplot,savename)
    if isnan(dS)
        %calculate expectation values
        v_ea_inc=expectation_val(length(timecell{1}),rho{1},sig{1,2});
        fig=figure,
        subplot(3,1,1)
        plot(timecell{1},fcell{1}.^2,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.2)
        xlabel('t');ylabel('f^2');
        subplot(3,1,2)
        plot(timecell{1},v_ea_inc.^2,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.2)
        xlabel('time(s)');ylabel('|<a>|^2');%title('<a>');
        subplot(3,1,3)
        plot(fcell{1}.^2,v_ea_inc.^2,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.2)
        xlabel('f^2');ylabel('|<a>|^2');%title('<a>');
        sgtitle(plottitle)
    else
        %calculate expectation values
        v_ea_inc=expectation_val(length(timecell{1}),rho{1},sig{1,2});
        v_ea_wait=expectation_val(length(timecell{2}),rho{2},sig{1,2});
        v_ea_dec=expectation_val(length(timecell{3}),rho{3},sig{1,2});
        %make plot
        fig=figure;
        subplot(3,1,1)
        plot(timecell{1},fcell{1}.^2,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.2)  %color dark red
        hold on
        plot(timecell{2},fcell{2}.^2,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.2)  %color green
        hold on
        plot(timecell{3},fcell{3}.^2,'color',[0, 0.4470, 0.7410],'LineWidth',1.2)       %color blue
        xlabel('t');ylabel('f^2');
        subplot(3,1,2)
        plot(timecell{1},v_ea_inc.^2,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.2)  %color dark red
        hold on
        plot(timecell{2},v_ea_wait.^2,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.2) %color green
        hold on
        plot(timecell{3},v_ea_dec.^2,'color',[0, 0.4470, 0.7410],'LineWidth',1.2)       %color blue
        xlabel('time(s)');ylabel('|<a>|^2');%title('<a>');
        subplot(3,1,3)
        plot(fcell{1}.^2,v_ea_inc.^2,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.2)  %color dark red
        hold on
        plot(fcell{2}.^2,v_ea_wait.^2,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.2) %color green
        hold on
        plot(fcell{3}.^2,v_ea_dec.^2,'color',[0, 0.4470, 0.7410],'LineWidth',1.2)       %color blue
        xlabel('f^2');ylabel('|<a>|^2');%title('<a>');
        sgtitle(plottitle)
    end
%     sgtitle(['JCH bistability: J=' num2str(J) ' \gamma=\kappa=' num2str(gamma) '  \mu=' num2str(mu) '  \epsilon=\omega - J'])
    if saveplot
        saveas(fig,savename)
    end

end