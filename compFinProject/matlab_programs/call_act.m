function y = call_act(S, T, K, r, sig, del, T_exp)
    for i=1:length(T)
        for j=1:length(S)
            if T(i) == T_exp
                y(i,j) = max(S(j) - K, 0);
            else
                d1 = (log(S(j)/K) + (r - del + 0.5*sig^2)*(T_exp - T(i)))/(sig*sqrt(T_exp - T(i)));
                d2 = d1 - sig*sqrt(T_exp - T(i));
                y(i,j) = S(j)*exp(-del*(T_exp-T(i))).*normcdf(d1) - K*exp(-r*(T_exp-T(i)))*normcdf(d2);
            end
        end
    end
end

