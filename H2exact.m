function OA = H2exact(f, clustconfig, sig )
%-------------------------------------------------------------------------%
%   Calculate the second part of the Hamiltonian depending on the time 
%   dependent coefficients, reducing it to single operators (matrices) 
%   of each cluster acting on their respective cluster
%   wave functions.
%Parameters:
%   f               pumping strength
%   clustconfig     matrix where the row index represents the cluster and
%                   contains the indices of the sites in that respective cluster.
%   NNM             Connection matrix, each row is a site on the lattice 
%                   and the columns with a 1 are the neirest neighbours 
%                   outside the cluster.
%   sig             contains the matrix operators for sigma x, y and z and
%                   a^dagger, a
%   SigExp          containt <a^dagger> and <a>
%-------------------------------------------------------------------------%
    n_sc = size(clustconfig, 2);            %number of sites in cluster
    n_coeff = size(sig{1,1},1);             %number of coefficients in cluster
    n_clust = size(clustconfig, 1);         %number of clusters

    OA = cell(1, n_clust);
    for i = 1:n_clust
        OA{i} = sparse(n_coeff,n_coeff);
%         OA{i}= OA{i} + f*(sig{1,1}+sig{1,2});

%         OA{i} = zeros(n_coeff,n_coeff);
    end

    for cluster=1:n_clust
        for i_clust = 1:n_sc                                                %index of the site inside the cluster.
                                                                            %this is just 1 loop clusters are size 1
%             i_site = clustconfig(cluster,i_clust);                          %index of the site in the lattice. 
            OA{cluster}= OA{cluster} + f*(sig{i_clust,1}+sig{i_clust,2});                                                                            
        end
        OA{cluster} = (-1i)*OA{cluster};        %factor -1i from schrodinger eq
    end
end