function iota = cap_prod_inv(opt, glob, p, psi)
    iota    = (exp((psi - p.delta) * p.kappa) - 1) / p.kappa + p.delta;
    
    % (iota - p.delta) .^ 2 - (iota - p.delta) + (psi - p.delta)     = 0;
    % iota = 1 / 2 - sqrt(1 - 4 * ((psi - p.delta))) / 2 + p.delta;
end