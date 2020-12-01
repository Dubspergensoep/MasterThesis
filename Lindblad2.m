function L=Lindblad2(A,Ad,Rho)
    %Ad=ctranspose(A);
    L=2*A*Rho*Ad-Ad*A*Rho-Rho*Ad*A;
end