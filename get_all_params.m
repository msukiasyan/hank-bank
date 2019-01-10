function pars = get_all_params(options, params)
    pars            = params;
    pars.K          = ((pars.r_minus + pars.delt) / pars.alp) ^ (1 / (pars.alp - 1));
    pars.w          = (1 - pars.alp) * pars.K ^ pars.alp;
    pars.alp_hat    = pars.f_bank / (pars.rho_bank + pars.f_bank - pars.r_F);
    pars.r_plus     = (pars.theta_bank * pars.r_F - pars.r_minus * pars.alp_hat) /...
                        (pars.theta_bank - pars.alp_hat);
    pars.x_a        = pars.alp_hat / pars.theta_bank;
    pars.lambda     = pars.alp_hat * (pars.r_minus - pars.r_plus) / pars.theta_bank;
end