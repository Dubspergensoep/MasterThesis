function [ time, Ct] = CalculateTrajectory(gamma, Kappa, J, f, w, epsilon, g, mu, Nxc, Nyc, Nx, Ny, T, dt, fill, Nmax, spin)
%-------------------------------------------------------------------------%
%Calculate the time evolution of the quantum trajectory.
%Parameters:
%   gamma       system parameter (spin dissipation)
%   J           hopping amplitude
%   f           pumping strength
%   Nxc         number of sites in x direction of cluster
%   Nyc         number of sites in y direction of cluster
%   Nx          number of rows
%   Ny          number of columns
%   T           Final time of quantum trajectory.
%   dt          Time steps to save.
%   wftype      type of wave function (C: Coherent, F: Fock)
%   fill        initial number of photons per site.
%   Nmax        Cutoff of local Hilbert space.
%   spin        spin vector for the generation of the initial wave function
%-------------------------------------------------------------------------%
    if (Nx ~= 1) && (Ny ~= 1)
        z = 4;                    %2d number of neirest neighbours
    else
        z = 2;                    %1d number of neirest neighbours
    end
    %Some system parameters and operators:
    N = Nx*Ny;                                                                                          %number of sites in lattice
    Nc = Nxc*Nyc;                                                                                       %number of sites in cluster
    NC = N/Nc;                                                                                          %number of clusters in lattice
    N_coeff = (2*(Nmax+1))^Nc;                                                                          %number of coefficients in (single) cluster wave function
     
    if Nxc == 1 && Nx ==1
        [~, ~, out] = CalcClusterNeighbours1D( Nxc, Nyc, Nx, Ny );
    else
        [~, ~, out] = CalcClusterNeighbours( Nxc, Nyc, Nx, Ny );                                        %calculates 'all' neirest neighbour as well as 'in' and 'out' of the cluster.
    end
    
    clustconfig = ClusterConfiguration(Nxc, Nyc, Nx, Ny);                                               %matrix where the row index represents the cluster and contains the indices of the sites in that respective cluster.
    sig = GetAllOperatorsSparse(Nc, Nmax);                                                              %calculate all the operators (sparse).
    NNM = NeirestNeighbourMatrix(out);                                                                  %Connection matrix, each row is a site on the lattice and the columns with a 1 are the neirest neighbours outside the cluster.
    A1 = H1(gamma, Kappa, J, 0, w, epsilon, g, mu, z, Nxc, Nyc, Nx, Ny, sig );                          %Time independt part of the Hamiltonian
    Tnorm = kron(eye(NC),ones(1,N_coeff));                                                              %Matrix to help calculate the norm of a product wave function.
    
    %Calculate the input wave function:
    Cin = GetFockWf(sig, N, fill, spin);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Save the wave function:
    C = transpose(Cin);
    
    %Continue time evolution of the trajectory
    for i = 1:floor(T/dt)
        disp(i)
        t_left = dt;                                                                                    %Time left untill next dt
        check = true;                                                                                   %Check variable to determine if timestep dt endend (Check = false).
        e1 = rand;                                                                                      %Random variable to determine jump time.
        ct=i*dt;                                                                                        %Current time needed for f(t)
        while check
            [t_s, Cin, ~, ~, ~] = EffH(J, f, z, Nxc, Nyc, Nx, Ny,...
                clustconfig, NNM, sig, A1, Cin, Tnorm, [0 t_left], e1,ct);                              %Time evolution up untill t_left or untill the jump time t_j.
            if t_s(end) == t_left                                                                       %NO JUMP in the chosen time interval.
                Cin = transpose(Cin(end,:))./NormC(transpose(Cin(end,:)), NC);                          %Normalize wave function
                C(size(C,1)+1,:) = transpose(Cin);                                                      %Save dt value.
                check = false;                                                                          %Check becomes false to go to the next time step.
            else                                                                                        %JUMP occurs on time t_j = t_s(end).
                i_jump = GetSite(ChanceInterval(transpose(Cin(end,:)), clustconfig, sig, Tnorm), rand); %Determine the site that jumps using their probability distribution.
                Cin = Jump(transpose(Cin(end,:)), clustconfig, sig, i_jump)./...
                    NormC(Jump( transpose(Cin(end,:)), clustconfig, sig, i_jump ), NC);                 %Update wave function with jump and normalize.
                t_left = t_left - t_s(end);                                                             %Time left untill next dt step.
                e1 = rand;                                                                              %Pick a new random number to determine jump time.
                ct=ct+t_s(end);                                                                         %Update current time for next loop
            end
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    time = (0:dt:floor(T/dt)*dt)';              %time interval
    Ct = C;                                     %wave function
    
end