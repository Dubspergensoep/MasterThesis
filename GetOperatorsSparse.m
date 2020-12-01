function A = GetOperatorsSparse( o , N)
%-------------------------------------------------------------------------%
%Calculates the operators on the product wave function for every site, 
%given the operator 'o' for that site.
%Parameters:
    %o      operator
    %N      number of sites in cluster
%-------------------------------------------------------------------------%
    o = sparse(o);
    
    %Size of the Hilbert space
    Hs = size(o,1);
    
    %data structures
    A = cell(1, N);
    %A = zeros(Hs^N,Hs^N,N);
    C = cell(1, N);
    
    %unity operators
    for i = 1:N
        C(i) = {speye(Hs)};
    end
    
    %apply kronecker product:
    for i = 1:N
        D = C;
        D(i) = {o};
        A{i} = superkronSparse(D);
    end
end