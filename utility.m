function u  = utility(c, opt, glob, p)
    u       = (max(c, 1e-12) .^ (1 - p.ga) - 1) / (1 - p.ga);
end