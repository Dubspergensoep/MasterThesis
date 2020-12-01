function [rho_lt,time]=Solve_drhodt(H1,sig,J,f,gamma,Kappa,Rho_0,dT,t_end)
    s_WF=length(Rho_0);
    %construct time vec
    t_e=ceil(t_end/dT)*dT; %make sure t_end is dividble by dT
    time=0:dT:t_e;
    
    %linearize rho_0
    Rho_l_0=Rho_0(:);
    
    % Constuct storage variables
    rho_lt=zeros(s_WF^2,length(time));
    rho_lt(:,1)=Rho_l_0;
    %tic
    for t=1:(length(time)-1)
        %construct time window
        tspan=[time(t) time(t+1)];
        %solve DV
        sol = ode15s(@(t,y) odefun_DM_MF(t, y,H1,sig,J, f, gamma, Kappa), tspan, Rho_l_0);
        %Updata Rho_l_0
        Rho_l_0=sol.y(:,end);
        %Save results
        rho_lt(:,t+1)=sol.y(:,end);
    end
    %toc
end