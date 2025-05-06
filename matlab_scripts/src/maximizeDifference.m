function [pairs, totalMaxDifference] = maximizeDifference(s, c)
    % Number of elements in each array
    n = length(s);

    % Create a list of edges for the graph, where each edge has (i, j, wt)
    edges = [];
    for i = 1:n
        for j = 1:n
            edges = [edges; i, j, c(i) - s(j)];  % Adjust for zero-based indexing
        end
    end
    
    % Call the maximum-weighted matching function
    % Assuming 'maxWeightMatching' is available and correctly implemented
    mate = maxWeightMatching(edges, true);

    % Extract pairs and calculate the total maximum difference
    pairs = zeros(n, 2);
    totalMaxDifference = 0;
    for i = 1:n
        if mate(i) ~= -1
            % Adjust for zero-based indexing in mate output
            pairs(i, :) = [i, mate(i)];  % +1 to convert back to MATLAB 1-based indexing
            totalMaxDifference = totalMaxDifference + (c(i) - s(mate(i)));
        end
    end
end
