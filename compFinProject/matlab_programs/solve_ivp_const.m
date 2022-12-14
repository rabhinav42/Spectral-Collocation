function [u,t, cpu_t] = solve_ivp_const(h0, opt, sig, r, E, T, m, N, D1, D2, s, ep, C)
    hn = h0;
    t = [0:hn:T];
    Xi = xi(N);
    S = [];
    for i=1:m
        s0 = s(i);
        s1 = s(i+1);
        S = [S X(Xi(1:N), s0, s1)];
    end
    S = [S s(end)];
    u = zeros(length(t), m*N+1);
    v = psi1(S, E, opt);
    u(1,:) = v;
    vn = v(2:end-1);
    v0 = vn;
    
    ern = 1;
    
    tic;
    for j=2:length(t)
        vn1 = vn + hn*F_call(t(j-1), vn, s, N, r, sig, D1, D2, E, m, opt, ep, C);
        if(opt == 'p')
            vn1 = max(vn1, v0);
        end
        vnm = vn + 0.5*hn*(F_call(t(j-1), vn, s, N, r, sig, D1, D2, E, m, opt, ep, C) + F_call(t(j), vn1, s, N, r, sig, D1, D2, E, m, opt, ep, C));
        if(opt == 'p')
            vnm = max(vnm, v0);
        end
        u(j,:) = [phi1(E, opt) vnm phi2(s(end), t(j), E, r, opt)];
        vn = vnm;
    end
    cpu_t = toc;
    
end
