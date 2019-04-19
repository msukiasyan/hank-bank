function pars = get_all_params(opt, p)
    p.K             = ((p.r_minus + p.delta) / p.Aprod / p.alpha) ^ (1 / (p.alpha - 1));
    p.w             = (1 - p.alpha) * p.Aprod * p.K ^ p.alpha;
    p.alp_hat       = p.f_bank / (p.rho_bank + p.f_bank - p.r_F);
    p.r_plus        = (p.theta_bank * p.r_F - p.r_minus * p.alp_hat) /...
                        (p.theta_bank - p.alp_hat);
    p.x_a           = p.alp_hat / p.theta_bank;
    p.lambda        = p.alp_hat * (p.r_minus - p.r_plus) / p.theta_bank;
    pars            = p;
end