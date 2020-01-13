function iota = cap_prod_inv(opt, glob, p, psi)
    iota    = (exp((psi - p.delta) * p.kappa) - 1) / p.kappa + p.delta;
end