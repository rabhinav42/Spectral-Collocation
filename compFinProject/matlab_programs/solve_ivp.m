function [u,t, cpu_t] = solve_ivp(tol, h0, opt, sig, r, E, T, m, N, D1, D2, s, ep, C)
    hn = h0;
    inttime = 0;
    t = [inttime];
    Xi = xi(N);
    S = [];
    for i=1:m
        s0 = s(i);
        s1 = s(i+1);
        S = [S X(Xi(1:N), s0, s1)];
    end
    S = [S s(end)];
    v = psi1(S, E, opt);
    u = v;
    vn = v(2:end-1);
    v0 = vn;
    
    ern = 1;
    
    tic;
    while(inttime < T)
        while(ern >= tol)
            vn1 = vn + hn*F_call(inttime, vn, s, N, r, sig, D1, D2, E, m, opt, ep, C);
            if(opt == 'p')
                vn1 = max(vn1, v0);
            end
            vnm = vn + 0.5*hn*(F_call(inttime, vn, s, N, r, sig, D1, D2, E, m, opt, ep, C) + F_call(inttime+hn, vn1, s, N, r, sig, D1, D2, E, m, opt, ep, C));
            if(opt == 'p')
                vnm = max(vnm, v0);
            end
            ern = norm(vnm - vn1);
            hn = hn/2;
        end
        
        hn = hn*2;
        hprev = hn;
        inttime = inttime + hn;
        vprev = vn;
        vn = vnm;
        u = [u;[phi1(E, opt) vn phi2(s(end), inttime, E, r, opt)]];
        t = [t; inttime];
        hn1 = hn*(0.5*tol/ern)^(1/3);
        hn = min(hn1,2);
        ern = 1;
    end
    
    inttime = inttime - hprev;
    t = t(1:end-1);
    hn = T - inttime;
    vn1 = vprev + hn*F_call(inttime, vprev, s, N, r, sig, D1, D2, E, m, opt, ep, C);
    if(opt == 'p')
        vn1 = max(vn1, v0);
    end
    vnm = vprev + 0.5*hn*(F_call(inttime, vprev, s, N, r, sig, D1, D2, E, m, opt, ep, C) + F_call(inttime+hn, vn1, s, N, r, sig, D1, D2, E, m, opt, ep, C));
    if(opt == 'p')
        vn1 = max(vn1, v0);
    end
    cpu_t = toc;
    inttime = T;
    u = u(1:end-1,:);
    u = [u; [phi1(E, opt) vnm phi2(s(end), inttime, E, r, opt)]];
    t = [t;inttime];
end

