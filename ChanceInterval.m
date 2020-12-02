function P = ChanceInterval( cin, clustconfig, sig, T )
%-------------------------------------------------------------------------%
%   Calculates the chance intervals of the sites to make a jump.
%Parameters:
%   cin             Input wavefunction in coefficient notation
%   clustconfig     Cluster configuration
%   sig             Jump Operators and corresponding dissipation rate.
%   T               Matrix to help calculate the norm of a product wave function.
%-------------------------------------------------------------------------%
    
    n = numel(clustconfig);
    J = zeros(n,1);
    J2 = J;
    for i = 1:n
        J(i) = Norm(Jump(cin, clustconfig, sig, i), T);
        J2(i) = sum(J(1:i));
    end
%     for i = 1:n
%         J(i) = sig{i,2}*Norm(Jump(cin, clustconfig, sig, i), T);
%         J2(i) = sum(J(1:i));
%     end
    P = J2/sum(J);
    
    
    %Give possible warning and errors:
    if (abs(P(end) - 1) > 1e-3)
        warning('Sum of the chances deviates more than 1e-3 from 1.')
    elseif isnan(P(end))
        error('Chance interval contains NaN data type.')
    end
end

