function psi = GetFockWf(sig, Nsites, fill, spin)
%-------------------------------------------------------------------------%
%   Calculate the Fock state (cluster) wave function.
%   functions still works without adaptation but all particles are
%   generated with spin up also fill must be smaller then
%   size(sig{1,1},1)/2 obviously
%Parameters:
%   sig         contains the operators on the system
%   Nsites      number of sites/modes in the system
%   fill        number of particles at each site
%   spin        vector containing spin of the start configuration.
%               IMPORTANT: spacing must be done by ; and not , or space
%-------------------------------------------------------------------------%
    %some system parameters:
    Ns = size(sig, 1);                      %number of  sites in cluster
    Nc = Nsites/Ns;                         %number of clusters in lattice
    fill = round(fill);                     %make sure it's a rounded number
    snos= size(sig{1,1},1)/2;                   %size no spin wave function

    %construct the single cluster wave function:
    wfc = zeros(snos,1);
    wfc(1) = 1;          %This was we start in spin down
    wfc=kron(spin,wfc);
    
    if fill > 0
        for i = 1:Ns
            for j = 1:fill
                wfc = sig{i,1}*wfc;
            end
        end
    end
    %normalize
    wfc = wfc/sqrt(wfc'*wfc);
    
    %construct cluster gutzwiller wave function:
    psi = repmat(wfc, Nc, 1);    
end

