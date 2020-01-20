function q = q_from_iota(opt, glob, p, iota)
    q       = 1 + p.kappa * (iota - p.delta);
    % q       = 1 ./ (1 - 2 *(iota - p.delta));
end