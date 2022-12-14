function y = put_act(S, T, E, r, sig, steps)
    % get american put price for given values of S at t=0 maturity time T 
    % with the binomial model.
    for j=1:length(S)
        [~, op] = binprice(S(j), E, r, T, T/steps, sig, 0, 0, 0, 0);
        y(j) = op(1,1);
    end
end

