  
function u  = utility_prime(c, h, opt, glob, p)
        u       = (max(c -  p.disutil * (h .^ (1 + 1/p.varphi)) / (1+1/p.varphi), 1e-12) .^ ( - p.ga));
    
end