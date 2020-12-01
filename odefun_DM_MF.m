function dydt=odefun_DM_MF(t,rho_l,H1,sig,J, f, gamma, kappa)
    % rho_l             rho linear
    
    % construct rho
    s_WF=sqrt(length(rho_l));
    rho=reshape(rho_l, s_WF, s_WF);
    
    % construct hamiltonian
    H=H1+H_MF_simple2(t,rho,sig,f,J);
    
    % equation
    drho= -1i*Com(H,rho) + (gamma/2)*Lindblad2(sig{2},sig{1},rho) + (kappa/2)*Lindblad2(sig{6},sig{5},rho);
    
    % linearize drho
    dydt=drho(:);
end