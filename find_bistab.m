function [rho,time,v_f,dS]=find_bistab(H1,sig,J,gamma,Kappa,Rho_0,dT,rate,fmin,fmax,twait,resstep,extra_steps,method)
    %Function will look for a resonance while increasing f.
    %when the resonance is found the system will be able to relax.
    %After this relaxation f will be decreased again.
    %if no resonance is found and f max is reach the program stops and f
    %will not be decreased
    rho=cell(1,3);
    time=cell(1,3);
    v_f=cell(1,3);
    s_WF=length(Rho_0);
    %resstep=ceil(resstep/dT)*dT;
    %% increase f
    t_fmax=fmax/rate;
    t_fM=ceil(t_fmax/dT)*dT; %make sure t_end is dividble by dT (Solve_drhodt also does this)
    %construct time vec
    t_fmin=fmin/rate;
    t_fm=ceil(t_fmin/dT)*dT; %make sure t_end is dividble by dT (Solve_drhodt also does this)
    f=@(t) laser(t,rate,t_fM+dT*extra_steps,0,0);    %This could be done simpler
    [rho_inc,time_inc]=Solve_drhodt(H1,sig,J,f,gamma,Kappa,Rho_0,dT,t_fm);
    
    % Looking for resonance
    %while loop parameters
    noRes=true;
    tc=time_inc(end); %current time
    
    while noRes && tc+dT*resstep <= t_fM
        % CHECK FOR RESONANCE
        %calculate v_ea
        v_ea=expectation_val(length(time_inc),rho_inc,sig{1,2});
        %look for resonance
        if checkRes(v_ea,method,dT)
            noRes=false;
        end
        
        % CALCULATE NEXT SEGMENT
        f_temp=@(t) f(tc+t);
        Rho_0=reshape(rho_inc(:,end), s_WF, s_WF);
        [rho_temp,time_temp]=Solve_drhodt(H1,sig,J,f_temp,gamma,Kappa,Rho_0,dT,dT*resstep);
        rho_inc=[rho_inc rho_temp];
        time_inc=[time_inc time_inc(end)+time_temp];
        % UPDATE tc
        tc=tc+dT*resstep;
        if tc+dT*resstep > t_fM
            disp('fmax was reached')
        end
    end
    
    
    if not(noRes) %this only need to happen if a resonance f was found
        % EXTRA INCREASE
        % f should be increased a little further after the resonance sinse we
        % to got slighty past the resonance frequency and not stop in the
        % middle of our resonance. This could be automate by looking at the
        % derivative again but i will not do this for now
        f_temp=@(t) f(tc+t);
        Rho_0=reshape(rho_inc(:,end), s_WF, s_WF);
        [rho_temp,time_temp]=Solve_drhodt(H1,sig,J,f_temp,gamma,Kappa,Rho_0,dT,dT*extra_steps);
        rho_inc=[rho_inc rho_temp];
        time_inc=[time_inc time_inc(end)+time_temp];
        % UPDATE tc
        tc=tc+dT*extra_steps;
        % SAVE ALL F VALUES
        f_inc=time_inc;
        for i=1:length(time_inc)
            f_inc(i)=f(time_inc(i));
        end
        t_fmax_inc=tc;  %time at maximal f value while increasing
        % STORE DATA
        rho{1}=rho_inc;
        time{1}=time_inc;
        v_f{1}=f_inc;
        %% WAIT
        tw=ceil(twait/dT)*dT;
        %construct new Rho_0
        Rho_0=reshape(rho_inc(:,end), s_WF, s_WF);
        %construct f wait
        f=@(t) t_fmax_inc*rate;
        %solve drhodt
        [rho_wait,time_wait]=Solve_drhodt(H1,sig,J,f,gamma,Kappa,Rho_0,dT,tw);
        f_wait=time_wait;
        for i=1:length(time_wait)
            f_wait(i)=f(time_wait(i));
        end
        time_wait=time_wait+time_inc(end);
        % STORE DATA
        rho{2}=rho_wait;
        time{2}=time_wait;
        v_f{2}=f_wait;
        %% DECREASE f
        te=ceil(t_fmax_inc/dT)*dT;
        %construct new Rho_0
        Rho_0=reshape(rho_wait(:,end), s_WF, s_WF);
        %construct f dec
        f=@(t) laser(t_fmax_inc+tw+t,rate,t_fmax_inc,tw,rate);
        %solve drhodt
        [rho_dec,time_dec]=Solve_drhodt(H1,sig,J,f,gamma,Kappa,Rho_0,dT,te);
        f_dec=time_dec;
        for i=1:length(time_dec)
            f_dec(i)=f(time_dec(i));
        end
        time_dec=time_dec+time_wait(end);
        % STORE DATA
        rho{3}=rho_dec;
        time{3}=time_dec;
        v_f{3}=f_dec;
        % CALCULATE BISTAB SURFACE dS
        v_ea_inc=expectation_val(length(time_inc),rho_inc,sig{1,2});
        v_ea_dec=expectation_val(length(time_dec),rho_dec,sig{1,2});
        flip_fdec=flip(f_dec);
        flip_veadec=flip(v_ea_dec);
        dS=trapz(flip_fdec.^2,flip_veadec.^2)-trapz(f_inc.^2,v_ea_inc.^2);
    else
        dS=NaN;
        % SAVE ALL F VALUES
        f_inc=time_inc;
        for i=1:length(time_inc)
            f_inc(i)=f(time_inc(i));
        end
        t_fmax_inc=tc;  %time at maximal f value while increasing
        % STORE DATA
        rho{1}=rho_inc;
        time{1}=time_inc;
        v_f{1}=f_inc;
        % STORE DATA
        rho{2}=[];
        time{2}=[];
        v_f{2}=[];
        % STORE DATA
        rho{3}=[];
        time{3}=[];
        v_f{3}=[];
    end
    
end

function res=checkRes(v_ea,method,dT)
    if strcmp(method.str,'ischange')
        %method to be implemnted
        res=false;
    elseif strcmp(method.str,'der_abs')
        %use derivative criteria aswell as an absolute thershold for v_ea
        %calculate derivatie of v_ea
        deadT=gradient(v_ea)/dT;
        if max(deadT)/mean(deadT)>method.threshold && v_ea(end)^2>method.abs_threshold
            res=true;
        else
            res=false;
        end
    else
        %calculate derivatie of v_ea
        deadT=gradient(v_ea)/dT;
        if max(deadT)/mean(deadT)>method.threshold
            res=true;
        else
            res=false;
        end
    end


end