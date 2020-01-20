function mpk = prod_K(opt, glob, p, cap, lab)
    mpk     = p.Aprod .* p.alpha .* (cap ./ lab) .^ (p.alpha - 1);
end