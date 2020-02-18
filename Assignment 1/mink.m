function Y=mink(X,k) % outputs a matrix with the first column containing the min value and the second the index for that value, the rows are in ascendning order with regards to value.
    min_val = min(X);
    %[n,m] = size(X);
    n = max(size(X));
    Y = [];
    for l=1:k
        for i=1:n
            if X(i) == min_val
                Y = [Y;X(i),i];
                X(i) = Inf;
                min_val = min(X);
                break
            end
        end
    end
    Y = Y(1:k,:);
end