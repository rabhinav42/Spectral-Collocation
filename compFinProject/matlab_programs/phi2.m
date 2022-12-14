function y = phi2(S, tau, E, r, opt)
    if(opt == 'c')
        y = S - E*exp(-r*tau);
    else
        y = 0;
    end
end