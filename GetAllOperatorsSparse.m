function B = GetAllOperatorsSparse( N, Nmax )
%-------------------------------------------------------------------------%
%Calculates the representations of a_i, a^+_i, a^+_ia_ifor the different 
%sites in the cluster. 
%Parameters:
%   N       number of sites in cluster
%   Nmax    maximum number of particles per site
%-------------------------------------------------------------------------%
    ap = kron(eye(2),GetCreationOperator(Nmax));            %creation operator
    am = kron(eye(2),GetAnnihilationOperator(Nmax));        %annihilation operator
    apm = ap*am;                                            %particle number operator 
    sz = kron([1 0; 0 -1], eye(Nmax+1));                    %pauli matrix sigma_z Nmax+1 omdat we tellen vanaf  |0>
    sp = kron([0 1; 0 0], eye(Nmax+1));                     %ladder operator sigma_plus
    sm = kron([0 0; 1 0], eye(Nmax+1));                     %ladder operator sigma_minus
    B = cell(N,6);
    
    %representations for the different sites:
    apc = GetOperatorsSparse(ap, N);
    amc = GetOperatorsSparse(am, N);
    apmc = GetOperatorsSparse(apm, N);
    szc = GetOperatorsSparse(sz, N);
    spc = GetOperatorsSparse(sp, N);
    smc = GetOperatorsSparse(sm, N);
    
    for i = 1:N
        B{i,1} = apc{i};
        B{i,2} = amc{i};
        B{i,3} = apmc{i};
        B{i,4} = szc{i};
        B{i,5} = spc{i};
        B{i,6} = smc{i};
    end
end

