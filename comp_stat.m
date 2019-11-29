function cs_res = comp_stat(opt, glob, p, par, low, high, npoints)
    cs_grid                 = linspace(low, high, npoints);
    cs_sol                  = cell(npoints, 1);
    cs_stats                = cell(npoints, 1);
    rates_old               = p.r_init;
    
    parfor piter = 1:npoints
        p_copy                              = p;
        p_copy.(par)                        = cs_grid(piter);
        if strcmp(par, 'zfactor')
            p_copy                          = update_idio_risk(opt, glob, p_copy);
        end
        [cs_sol{piter}, cs_stats{piter}]    = find_ss(opt, glob, p_copy, rates_old);
    end
    cs_res.cs_sol           = cs_sol;
    cs_res.cs_stats         = cs_stats;
    cs_res.cs_grid          = cs_grid;
    cs_res.par              = par;
end