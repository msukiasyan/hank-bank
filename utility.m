function u  = utility(c, opt, glob, p)
    u       = (c .^ (1 - p.ga) - 1) / (1 - p.ga);
end