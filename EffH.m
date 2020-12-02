function [time, Ct, tj, Cj, ij] = EffH( J, f, z, Nxc, Nyc, Nx, Ny, clustconfig, NNM, sig, A1, Cin, T, Tspan, eps, ct )
%-------------------------------------------------------------------------%
%Function EffH 
%   Calculates the time evolution of the wave function.
%Parameters:
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
%   sig             contains the matrix operators for sigma a, a^+, a^+a.
%   A1              matrix containing the Hamiltonian independent of
%                   time, per cluster.
%   Cin             input wave function.
%   T               matrix which makes the sum over the moduli of the
%                   coefficents of the different cluster wave functions
%   Tspan           Time evolution up till time T
%   eps             the random value to determine the jump time
%   ct              current time, this is need for f(t)
%-------------------------------------------------------------------------%
options = odeset('events', @(t,C) myEvent(t, C, eps, T));
[time, Ct, tj, Cj, ij] = ode45(@(t, C) odeC(t, C, J, f, z, Nxc, Nyc, Nx, Ny, clustconfig, NNM, sig, A1, ct), Tspan, Cin, options);
end