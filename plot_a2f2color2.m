function plot_a2f2color2(v_ea,time,plottitle,saveplot,indices)
    fig=figure,
    plot(time(1:indices(1)),v_ea(1:indices(1)).^2,'color','red')
    hold on
    plot(time((indices(1)+1):indices(2)),v_ea((indices(1)+1):indices(2)).^2,'color','green')
    hold on
    plot(time((indices(2)+1):end),v_ea((indices(2)+1):end).^2,'color','blue')
    xlabel('time step');ylabel('|<a>|^2');%title('<a>');
    legend('increasing f','constant f', 'decreasing f')
    title(plottitle)
%     sgtitle(['JCH bistability: J=' num2str(J) ' \gamma=\kappa=' num2str(gamma) '  \mu=' num2str(mu) '  \epsilon=\omega - J'])
    if saveplot.saveplot
        saveas(fig,saveplot.savename)
    end

end