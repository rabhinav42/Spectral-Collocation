function [R1, R2] = R(x, r, sig, s0, s1)
    % evaluate R1(x) and R2(x)
    y = X(x, s0, s1);
    R1 = r*y*2./(s1 - s0);
    R2 = 0.5*(sig^2)*(2*y./(s1 - s0))^2;
end

