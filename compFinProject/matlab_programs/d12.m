function [d1, d2] = d12(s,r,sig,tau)
    d1 = (log(s) + (r + 0.5*sig^2)*tau)/(sig*sqrt(tau));
    d2 = d1 - sig*sqrt(tau); 
end

