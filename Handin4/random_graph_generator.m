function W = random_graph_generator(n,k)
% Generating a random graph according to the preferential attachment model

k0 = k+1; % Number of nodes in current graph, k0 > c
W = ones(k0,k0)-diag(ones(k0,1)); % Adjacency matrix
W = sparse(W); % Transform W into a sparse matrix

for i=(k0+1):n
    if mod(i,2)==0 % If i even number
        c = floor(k/2); % Degree of new node
    else % If i odd number
        c = ceil(k/2); % Degree of new node
    end
    w = sum(W,2); % Degree vector of graph
    P = w./sum(w); % Probability vector for adding links
    for j=1:c % Choose c neighbours
        % Choose one neighbour with prob. prop. to node degrees
        neighbour = randsample(1:k0,1,true,full(P));
        % Assure not to choose the same neighbour again
        P(neighbour) = 0;
        W(k0+1,neighbour) = 1; % Add link (one direction)
        W(neighbour,k0+1) = 1; % Add link (other direction)
    end
    k0 = k0 + 1; % Number of nodes in current graph
end

end