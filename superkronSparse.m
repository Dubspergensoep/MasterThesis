function M = superkronSparse(D)
    % D is a cell array containing the operators and the order in which
    % the tensor product should be taken (we take the product from right to left).
    
    if length(D) == 1
        M = D{1};        
    else
        M = kron(D{end - 1}, D{end}); 

        if length(D) > 2
            for i=2:length(D)-1
                M = kron(D{end-i}, M);
            end
        end
    end
end