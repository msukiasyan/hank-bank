function u  = utility_prime(c, opt, glob, p)
    u       = max(c, 1e-12) .^ (- p.ga);
end