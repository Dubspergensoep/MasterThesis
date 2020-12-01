function a = GetAnnihilationOperator( Nmax )
%-------------------------------------------------------------------------%
%Calculates the representations of the annihilation operator for a system where
%the number of particles has a maximum Nmax for a certain mode.
%Note: the wave function/column vector is presumed to have been indexed 
%      from a state |0> to a state |Nmax>, where |0> is the first('highest') 
%      entry of the column vector and |Nmax> the last('lowest').
%Parameters:
%   Nmax   maximum number of particles in mode.
%-------------------------------------------------------------------------%
    %defining the dimensions of the operator
    Ncoeff = Nmax + 1;
    a = zeros(Ncoeff);
    
    %construct the operator
    for i = 1:Ncoeff-1
        a(i,i+1) = sqrt(i);        
    end
end