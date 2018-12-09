function pars = get_all_params(r_F, r_minus)
    global Aprod alp delt rho_bank f_bank theta_bank;
    pars.r_F        = r_F;
    pars.r_minus    = r_minus;
    pars.K          = ((r_minus + delt) / alp) ^ (1 / (alp - 1));
    pars.alp_hat    = f_bank / (rho_bank + f_bank - r_F);
    pars.r_plus     = (theta_bank * r_F - r_minus * pars.alp_hat) / (theta_bank - pars.alp_hat);
    pars.x_a        = pars.alp_hat / theta_bank;
    pars.lambda     = pars.alp_hat * (r_minus - pars.r_plus) / theta_bank;
end