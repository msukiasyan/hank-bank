function psi = cap_prod(opt, glob, p, iota)
    psi     = 1 / p.kappa * log(1 + p.kappa * (iota - p.delta)) + p.delta;
    % psi     = (iota - p.delta) - (iota - p.delta) .^ 2 + p.delta;
end