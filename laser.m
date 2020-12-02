function f=laser(t,r_i,t_fmax,tw,r_d)
    if t<=t_fmax
        f=r_i*t;
    elseif t<=t_fmax+tw
        f=r_i*t_fmax;
    elseif t>t_fmax+tw
        f=r_i*t_fmax-r_d*(t-(t_fmax+tw));
        if f<0.0
            f=0.0;
        end
    end
end