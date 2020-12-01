function plot_a2f2(time,v_ea,f_v,plottitle,saveplot)
    fig=figure,
    subplot(3,1,1)
    plot(time,f_v.^2,'color','black')
    xlabel('t');ylabel('f^2');
    subplot(3,1,2)
    plot(time,v_ea.^2,'color','black')
    xlabel('time(s)');ylabel('|<a>|^2');%title('<a>');
    subplot(3,1,3)
    plot(f_v.^2,v_ea.^2,'color','black')
    xlabel('f^2');ylabel('|<a>|^2');%title('<a>');
    sgtitle(plottitle)
%     sgtitle(['JCH bistability: J=' num2str(J) ' \gamma=\kappa=' num2str(gamma) '  \mu=' num2str(mu) '  \epsilon=\omega - J'])
    if saveplot.saveplot
        saveas(fig,saveplot.savename)
    end

end