function y = F_call(tau, u, s, N, r, sig, D1, D2, E, m, opt, ep, C)
    
% assuming u does not have the boundary terms, ie is of size mN-1
    
    y = [];
    % F1,j
    i = 1;
    Xi = xi(N);
    Zi = Z(s(1), s(2), N, r, sig, D1, D2);
    ynew = Zi(2:N+1, 2:N+1)*u(1:N)' + Zi(2:N+1,1)*phi1(E,opt) + g(u(1:N), E, X(Xi(2:N+1), s(i), s(i+1)), ep, C, opt)';
    %ynew = Zi(2:N+1, 2:N+1)*u(1:N)';
    y = [y ynew'];

    
    
    %Fi,j for i=2,...,m-1
    du1 = du(D1, [phi1(E,opt) u(1:N-1)], u(N+1:2*N), N);
    %du1 = du(D1, [0 u(1:N-1)], u(N+1:2*N), N);
    for i=2:m-1
        Zi = Z(s(i), s(i+1), N, r, sig, D1, D2);
        du0 = du1;
        if i==m-1
            du1 = du(D1, u(((i-1)*N): (i*N -1)), [u(((i)*N + 1) : ((i+1)*N-1)) phi2(s(end), tau, E, r, opt)], N);
            %du1 = du(D1, u(((i-1)*N): (i*N -1)), [u(((i)*N + 1) : ((i+1)*N-1)) phi_c(s(end), tau, E, r)], N);
        else
            du1 = du(D1, u(((i-1)*N): (i*N - 1)), u(((i)*N + 1) : ((i+1)*N)), N);
        end
        ynew = Zi(2:N+1, 1)*du0 + Zi(2:N+1, 2:N)*u(((i-1)*N + 1) : (i*N - 1))' + Zi(2:N+1, N+1)*du1 + g(u(((i-1)*N + 1) : i*N), E, X(Xi(2:N+1), s(i), s(i+1)), ep, C, opt)';
        %ynew = Zi(2:N+1, 1)*du0 + Zi(2:N+1, 2:N)*u(((i-1)*N + 1) : (i*N - 1))' + Zi(2:N+1, N+1)*du1;
        y = [y ynew'];
    end
        
    
    %Fm,j
    i = m;
    Zi = Z(s(i), s(i+1), N, r, sig, D1, D2);
    du0 = du1;
    ynew = Zi(2:N, 1)*du0 + Zi(2:N, 2:N)*u(((i-1)*N + 1) : (i*N - 1))' + Zi(2:N, N+1)*phi2(s(end), tau, E, r, opt) + g(u(((i-1)*N + 1) : (i*N - 1)), E, X(Xi(2:N), s(i), s(i+1)), ep, C, opt)';
    %ynew = Zi(2:N, 1)*du0 + Zi(2:N, 2:N)*u(((i-1)*N + 1) : (i*N - 1))' + Zi(2:N, N+1)*phi_c(s(end), tau, E, r);
    y = [y ynew'];
end