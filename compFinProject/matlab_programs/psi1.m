function y = psi1(S,E,opt)
    if(opt == 'c' || opt == 'b')
        y = max(S-E,0);
    else
        y = max(E-S,0);
    end
end