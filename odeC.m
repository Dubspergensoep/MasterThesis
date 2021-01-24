function dC = odeC(t, C, J, f, z, Nxc, Nyc, Nx, Ny, clustconfig, NNM, sig, A1, ct)
%-------------------------------------------------------------------------%
%   system of differential equations for the coefficients of the individual
%   site wave functions.
%Parameters:
%   t               time variable
%   C               coefficient value the ODE will calculate (column vector).
%   gamma           system parameter
%   J               hopping amplitude
%   f               pumping strength
%   z               number of neighbours
%   Nxc             number of sites in x direction of cluster
%   Nyc             number of sites in y direction of cluster
%   Nx              number of rows
%   Ny              number of columns
%   clustconfig     matrix where the row index represents the cluster and
%                   contains the indices of the sites in that respective cluster.
%   NNM             Connection matrix, each row is a site on the lattice 
%                   and the columns with a 1 are the neirest neighbours 
%                   outside the cluster.
%   sig             contains the matrix operators for sigma x, y and z.
%   switch_xy,
%   magn_y, magn_z  Alternative way to calculate expectation values, uses
%                   the permutation properties of the pauli matrices.
%   A1              matrix containing the Hamiltonian independt of
%                   time, per cluster.
%   ct              current time, this is need for f(t)
%-------------------------------------------------------------------------%
    Nc = Nxc*Nyc;                               %number of sites in cluster.
    N=Nx*Ny/(Nc);                               %number of clusters in lattice.  
    Ncoeff = size(sig{1,1},1);                  %number of coefficients in a cluster wave function.
    
    
    if (Nxc == Nx) && (Nyc == Ny)%exact solution
%         A2 = cell(1);
        A2 = H2exact(f(ct+t),clustconfig, sig);
    else%cluster solution
        SigExp = CalcExpValues( C, clustconfig, sig );       %We calculate the expectation value of sigma_i^r, with r= x, y, z.
        %f_t=f(ct+t);                                        
        A2 = H2( J, f(ct+t), z, clustconfig, NNM, sig, SigExp );             %Time dependent part of the Hamiltonian.
    end
                                                                            
    dC = zeros(N*Ncoeff,1);
    for j=1:N      %sum over the clusters.
        range = (Ncoeff)*(j-1)+1:(Ncoeff)*(j-1)+1+(Ncoeff-1);
        if (Nxc == Nx) && (Nyc == Ny)
            dC(range) = (A1{j}+A2{j})*C(range);
        else
            dC(range) = (A1{j}+A2{j})*C(range);
        end
    end
end   

