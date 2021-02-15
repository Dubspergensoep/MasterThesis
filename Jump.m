function c = Jump( cin, clustconfig, c_ops, i )
%-------------------------------------------------------------------------%
%   Gives the wave function when a jump occured on site i
%Parameters:
%   cin             wave function in coefficient notation.
%   sig             Cluster wave function from cluster containing site i,
%                   containing the spin up, spin down states.
%   clustconfig     configuration of the lattice sites in the cluster.
%   i               index of site where the jump occurs.
%-------------------------------------------------------------------------%
n_coeff = size(c_ops{1,1},1);                                           %number of coefficients in cluster wave function.
i_site = Mod(i, numel(clustconfig),1);                                  %Adjust i to to correct value
[x, y] = find(clustconfig == i_site);                                   %x = cluster index and y = site index in cluster.
range = n_coeff*(x-1)+1:n_coeff*(x-1)+1+(n_coeff-1);                    %Denotes the indices of the cluster wave function in the system wave function.
c_ops_i = y + floor((i-1)/numel(clustconfig))*round(length(c_ops)/2);   %Select lowering operator
cin(range) = c_ops{c_ops_i}*cin(range);                                 %Lowering operator on site i works on the wave function of the respective cluster.
c = cin;

end

