function y = X(x, s0, s1)
    % transform [-1,1] to [s0, s1]
    % i.e., X(x, s0, s1) = S where S is in [s0, s1] and x in [-1,1]
    y = 0.5*(s1 - s0).*x + 0.5*(s1 + s0);
end

