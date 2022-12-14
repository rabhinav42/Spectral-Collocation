function y = uno_act(S, T, K, r, sig, B, T_exp)
    % assuming del = 0 and S <= B
    for i=1:length(T)
        tau = T_exp - T(i);
        for j=1:length(S)
            if(T(i) == T_exp)
                y(i,j) = max(S(j) - K,0);
            else
                [d1B, d2B] = d12(S(j)/B, r, sig, tau);
                [d1BB, d2BB] = d12(B/S(j), r, sig, tau);
                [d1K, d2K] = d12(S(j)/K, r, sig, tau);
                [d1BK, d2BK] = d12(B*B/(K*S(j)), r, sig, tau);
                y(i,j) = S(j)*(normcdf(d1K) - normcdf(d1B)) - ...
                    K*exp(-r*tau)*(normcdf(d2K) - normcdf(d2B)) - ...
                    B*(S(j)/B)^(-2*r/sig^2)*(normcdf(d1BK) - normcdf(d1BB)) + ...
                    K*exp(-r*tau)*(S(j)/B)^(-2*r/sig^2 + 1)*(normcdf(d2BK) - normcdf(d2BB));
            end
        end
end