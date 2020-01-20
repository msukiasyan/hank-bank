function u  = utility(c, h, opt, glob, p)
    %u       = (max(c, 1e-12) .^ (1 - p.ga) - 1) / (1 - p.ga);
    u       = (max(c - p.disutil * (h .^ (1 + 1/p.frisch)) / (1+1/p.frisch), 1e-12) .^ (1 - p.ga) - 1) / (1 - p.ga);
    
end