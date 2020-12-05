function [rho,timecell,fcell,dS]=SSP(mu,g,J,N_max,dT,rate,fmin,fmax,twait,resstep,extra_steps,method,vplot,saveloc,saveplot,savelocplots,odemethod,dispt)
    %Solve Save and Plot
    gamma=0.005;
    Kappa=0.005;
    %J=0.126;
    omega=1;
    eps=omega-J;
    
    % sym parameters
    fill=0;
    spin=[0;1];
%     N_max=5;
    s_WF=2*(N_max+1);       %size wave function
    sig = GetAllOperatorsSparse(1,N_max);
    
    %Construct rho_0
    Psi=GetFockWf(sig,1,fill,spin);
    Rho_0=Psi*Psi';
    Rho_l_0=Rho_0(:);
    
    % Construct H1
    H1=H_MF_simple1(sig,omega,mu,eps,g);
    
    % Savename
    filename=[saveloc 'SSP_mu' num2str(mu) '_g' num2str(g) '_J' num2str(J)...
        '_fmin' num2str(fmin) '_fmax'  num2str(fmax) '_rate' ...
        num2str(rate) '_twait' num2str(twait) '_' num2str(dT) '_resstep' num2str(resstep) ...
        '_extra_steps' num2str(extra_steps) '_method' method.str ...
        '_' num2str(method.threshold) '_' num2str(method.abs_threshold) '.mat'];
    %%
    dS=zeros(1,2);  %dS(1) = surface in <n(f)> plot
                    %dS(2) = surface in <a(f)> plot
    %% Solve ODE
    if isfile(filename)
        load(filename)
    else
        [rho,timecell,fcell,dS(2)]=find_bistab(H1,sig,J,gamma,Kappa,Rho_0,dT,rate,fmin,fmax,twait,resstep,extra_steps,method,true,odemethod,dispt);
        save(filename,'rho','timecell','fcell','dS');
    end
    %% Calculate dSn
    % CALCULATE BISTAB SURFACE dS
    v_n_inc=expectation_val(length(timecell{1}),rho{1},sig{1,3});
    v_n_dec=expectation_val(length(timecell{3}),rho{3},sig{1,3});
    flip_fdec=flip(fcell{3});
    flip_vndec=flip(v_n_dec);
    dS(1)=trapz(flip_fdec.^2,flip_vndec.^2)-trapz(fcell{1}.^2,v_n_inc.^2);
    %% Plots
    savenameplot=[savelocplots 'SSP_plot_mu' num2str(mu) '_g' num2str(g) '_J' num2str(J)...
        '_fmin' num2str(fmin) '_fmax'  num2str(fmax) '_rate' ...
        num2str(rate) '_twait' num2str(twait) '_' num2str(dT) '_resstep' num2str(resstep) ...
        '_extra_steps' num2str(extra_steps) '_method' method.str ...
        '_' num2str(method.threshold) '_' num2str(method.abs_threshold) '.png'];
    plottitle=['JCH bistability: J=' num2str(J) '  \mu=' num2str(mu) 'g=' num2str(g) '\epsilon=\omega - J'];
    if length(vplot)==3
        %more plot option can be added later like in SaveandPlot
        if vplot(1)
            parforplot_a2f2color(rho,timecell,fcell,dS,sig,plottitle,saveplot,savenameplot);
        end
        if vplot(2)
        %calculate expectation values
        v_n_inc=expectation_val(length(timecell{1}),rho{1},sig{1,3});
        v_n_wait=expectation_val(length(timecell{2}),rho{2},sig{1,3});
        v_n_dec=expectation_val(length(timecell{3}),rho{3},sig{1,3});
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
        plot(timecell{1},v_n_inc.^2,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.2)  %color dark red
        hold on
        plot(timecell{2},v_n_wait.^2,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.2) %color green
        hold on
        plot(timecell{3},v_n_dec.^2,'color',[0, 0.4470, 0.7410],'LineWidth',1.2)       %color blue
        xlabel('time(s)');ylabel('|<a*a>|');%title('<a>');
        subplot(3,1,3)
        plot(fcell{1}.^2,v_n_inc.^2,'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.2)  %color dark red
        hold on
        plot(fcell{2}.^2,v_n_wait.^2,'color',[0.4660, 0.6740, 0.1880],'LineWidth',1.2) %color green
        hold on
        plot(fcell{3}.^2,v_n_dec.^2,'color',[0, 0.4470, 0.7410],'LineWidth',1.2)       %color blue
        xlabel('f^2');ylabel('|<a*a>|');%title('<a>');
        sgtitle(plottitle)
        end
    else
%         disp(class(timecell))
%         disp(class(rho))
        parforplot_a2f2color(rho,timecell,fcell,dS,sig,plottitle,saveplot,savenameplot);
    end
end