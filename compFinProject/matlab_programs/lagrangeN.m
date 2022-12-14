function y = lagrangeN(N, j, x)
    % evaluate L_j at x where we assume j is from 1, ..., N+1
    Xi = xi(N);
    z = Xi(Xi ~= Xi(j));
    denom = prod(Xi(j) - z);
    num = prod(x - z);
    y = num/denom;
end

