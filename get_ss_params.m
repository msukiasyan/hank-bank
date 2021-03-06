function pars = get_ss_params(opt, glob, p)
    p               = update_idio_risk(opt, glob, p);
    p.N             = sum(p.z .* p.z_dist);
    p.K             = ((p.r_minus + p.delta) / (p.Aprod * p.alpha)) ^ (1 / (p.alpha - 1)) * p.N;
    p.w             = (1 - p.alpha) * p.Aprod * p.K ^ p.alpha * p.N ^ (- p.alpha);
    p.eta           = p.f_bank / (p.rho_bank + p.f_bank - p.r_F);
    p.r_plus        = (p.theta_bank * p.r_F - p.r_minus * p.eta) /...
                        (p.theta_bank - p.eta);
    p.x_a           = p.eta / p.theta_bank;
    p.lambda        = p.eta * (p.r_minus - p.r_plus) / p.theta_bank;
    p.spread        = p.r_minus - p.r_plus;
    
    if opt.GK
        p.NW        = p.MGK / (p.f_bank - p.r_F);
    end
    
    pars            = p;
end