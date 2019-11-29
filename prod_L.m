function mpk = prod_L(opt, glob, p, cap)
    mpk     = p.Aprod .* (1 - p.alpha) .* cap .^ (p.alpha);
end