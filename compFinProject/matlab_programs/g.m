function y = g(v, E, S, ep, C, opt)
    if(opt == 'p')
        y = ep*C./(v + ep - E + S);
    else
        y = 0;
    end
end