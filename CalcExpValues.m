function SigExp = CalcExpValues( C, clustconfig, sig)
%-------------------------------------------------------------------------%
%   Calculate the expectation value of the x,y,z magnetization of each
%   site.
%Parameters:
%   C               Coefficient value the ODE will calculate
%   clustconfig     matrix where the row index represents the cluster and
%                   contains the indices of the sites in that respective 
%                   cluster.
%   sig             operators
%-------------------------------------------------------------------------%
    n_sc = size(clustconfig, 2);            %number of sites in cluster
    n_coeff = size(sig{1,1},1);             %number of coefficients in cluster
    n_clust = size(clustconfig, 1);         %number of clusters
%     sigexp = zeros(n_clust*n_sc, 1);      %This caused errors 1 should be
                                            %2 right?
    sigexp = zeros(n_clust*n_sc, 2);
    
    for i = 1:n_clust                       %Sum over clusters
        range = n_coeff*(i-1)+1:n_coeff*(i-1)+1+(n_coeff-1);
        d = C(range);
        for j = 1:n_sc                      %Sum over sites in cluster
            Ind = clustconfig(i,j);
            
            sigexp(Ind, 1) = d'*sig{j,1}*d;                 %d'*a^+*d
            
            sigexp(Ind, 2) = d'*sig{j,2}*d;                 %d'*a*d
        end
    end
    SigExp = sigexp;
end