function q = q_from_iota(opt, glob, p, iota)
    q       = 1 + p.kappa * (iota - p.delta);
end