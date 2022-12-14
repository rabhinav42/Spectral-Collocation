function y = xi(N)
    % return CGL points
    k = [0:N];
    y = flip(cos(k*pi/N));
end

