function [rho,timecell,fcell,dS]=SSP(mu,g,J,dT,rate,fmin,fmax,twait,resstep,extra_steps,method,vplot,saveloc,saveplot,savelocplots)
    %Solve Save and Plot
    gamma=0.005;
    Kappa=0.005;
    %J=0.126;
    omega=1;
    eps=omega-J;
    
    % sym parameters
    fill=0;
    spin=[0;1];
    N_max=5;
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
    %% Solve ODE
    if isfile(filename)
        load(filename)
    else
        [rho,timecell,fcell,dS]=find_bistab(H1,sig,J,gamma,Kappa,Rho_0,dT,rate,fmin,fmax,twait,resstep,extra_steps,method);
        save(filename,'rho','timecell','fcell','dS');
    end
    %% Plots
    savenameplot=[savelocplots 'SSP_plot_mu' num2str(mu) '_g' num2str(g) '_J' num2str(J)...
        '_fmin' num2str(fmin) '_fmax'  num2str(fmax) '_rate' ...
        num2str(rate) '_twait' num2str(twait) '_' num2str(dT) '_resstep' num2str(resstep) ...
        '_extra_steps' num2str(extra_steps) '_method' method.str ...
        '_' num2str(method.threshold) '_' num2str(method.abs_threshold) '.png'];
    plottitle=['JCH bistability: J=' num2str(J) '  \mu=' num2str(mu) 'g=' num2str(g) '\epsilon=\omega - J'];
    if length(vplot)==3
        %more plot option can be added later like in SaveandPlot
        parforplot_a2f2color(rho,timecell,fcell,dS,sig,plottitle,saveplot,savenameplot);
    else
%         disp(class(timecell))
%         disp(class(rho))
        parforplot_a2f2color(rho,timecell,fcell,dS,sig,plottitle,saveplot,savenameplot);
    end
end