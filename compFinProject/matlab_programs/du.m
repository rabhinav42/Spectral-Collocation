function y = du(D1, u0, u1, N)
    % evaluate dui(tau)
    denom = D1(N+1, N+1) - D1(1,1);
    y = (u1*D1(2:N+1, 1) - u0*D1(1:N, N+1))/denom;
end

