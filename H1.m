function OA = H1(gamma, Kappa, J, f, w, epsilon, g, mu, z, Nxc, Nyc, Nx, Ny, sig )
%-------------------------------------------------------------------------%
%   Calculate the first part of the Hamiltonian, the time independent part,
%   reducing it to single operators (matrices) of each cluster acting on 
%   their respective cluster wave functions.
%Parameters:
%   gamma           system parameter
%   J               hopping amplitude
%   f               pumping strength (This is moved to H2)
%   w               energy bosonic mode
%   epsilon         transition energy 2 level system
%   g               coupling photon and 2 level system
%   mu              chemical potential
%   z               number of neighbours
%   Nxc             number of sites in x direction of cluster
%   Nyc             number of sites in y direction of cluster
%   Nx              number of rows
%   Ny              number of columns
%   sig             operators
%-------------------------------------------------------------------------%
    %Some system parameters:
    c = ClusterConfiguration(Nxc, Nyc, Nx, Ny);
    SiteI = SiteClusterIndex(Nxc, Nyc, Nx, Ny);
    if Nxc == 1 && Nx == 1
        [~, in, ~] = CalcClusterNeighbours1D(Nxc, Nyc, Nx, Ny);
    else
        [~, in, ~] = CalcClusterNeighbours(Nxc, Nyc, Nx, Ny);
    end

    %The calculation:
    OA = cell(1,size(c,1));
    for cluster=1:size(c,1)
        A = sparse(size(sig{1,1},1), size(sig{1,1},2));                     
%         A = zeros(size(sig{1,1}));
        for i_clust = 1:length(c(cluster,:))                                    %index of the site inside the cluster.
            i_site = c(cluster,i_clust);                                        %index of the site in the lattice.
            
            % geen 1i bijdrage hier doordat deze wegvalt bij het product
            % met die -1i uit de Schrodinger vergelijking.
            A = A + (-gamma/2)*sig{i_clust,3};
            A = A + (-Kappa/2)*sig{i_clust,5}*sig{i_clust,6};

            %Bijdrage van (-1i) komt van de schrodinger vergelijking.
            A = A + (-1i)*((w - mu)*sig{i_clust, 3} ...
                + 0.5*(epsilon-mu)*sig{i_clust,4} ...
                + g*(sig{i_clust,2}*sig{i_clust,5} + sig{i_clust,1}*sig{i_clust,6})); % Jaynes-Cummings part
            %A = A + (-1i)*f*(sig{i_clust,1}+sig{i_clust,2});                         % Term corresponding to photon pumping
            
            n_in = nonzeros(in(i_site,:));      %neighbours of i_site in cluster. !in lattice picture!

            for j = 1:length(n_in)
                ind = SiteI(n_in(j));                                           %Cluster index of the site n_in(j);
                
                %Hopping contribution inside cluster:
                A = A + (1i)*(J/z)*(sig{i_clust,1}*sig{ind, 2});
               
            end
        end
        OA{cluster} = sparse(A);
    end
end