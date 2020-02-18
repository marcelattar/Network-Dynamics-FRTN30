function Y=maxk(X,k) % outputs a matrix with the first column containing the max value and the second the index for that value, the rows are in descending order with regards to value.
    max_val = max(X);
    %[n,m] = size(X');
    n = max(size(X));
    Y = [];
    for l=1:k
        for i=1:n
            if X(i) == max_val
                Y = [Y;X(i),i];
                X(i) = 0;
                max_val = max(X);
                break
            end
        end
    end
    Y = Y(1:k,:);
end