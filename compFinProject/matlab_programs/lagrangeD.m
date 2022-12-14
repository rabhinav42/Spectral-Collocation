function [D1, D2] = lagrangeD(N)
    % evaluate D1 and D2
    % that is first and second derivatives of L_j at xi_k for each j and k
    [D1,~] = cheb(N);
    D1 = D1';
    D2 = D1*D1;
end

