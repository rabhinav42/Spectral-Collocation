function Zi = Z(s0, s1, N, r, sig, D1, D2)
    % evaluate Z_(i,j,k)
    Xi = xi(N);
    for j=1:N+1
        [R1, R2] = R(Xi(j), r, sig, s0, s1);
        for k=1:N+1
            L = lagrangeN(N,k,Xi(j));
            Zi(j,k) = R1*D1(k,j) + R2*D2(k,j) - r*L;
        end
    end        
end

