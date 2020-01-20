function mpk = prod_L(opt, glob, p, cap, lab)
    mpk     = p.Aprod .* (1 - p.alpha) .* (cap ./ lab) .^ (p.alpha);
end