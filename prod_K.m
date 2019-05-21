function mpk = prod_K(opt, glob, p, cap)
    mpk     = p.Aprod * p.alpha * cap .^ (p.alpha - 1);
end