function H=H_MF_simple2(t,rho,sig,f,J)
    %ft=0.00001*t;
    %time dependent
    ea=trace(rho*sig{2});   %<a>
    ead=trace(rho*sig{1});  %<a*>
    H = f(t)*(sig{1}+sig{2})- J*(ea*sig{1}+ead*sig{2});
end