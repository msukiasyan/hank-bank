function cs_res = comp_stat(opt, p, par, low, high, npoints)
    cs_grid                 = linspace(low, high, npoints);
    cs_sol                  = cell(npoints, 1);
    cs_stats                = cell(npoints, 1);
    rates_old               = [0.038 0.04];
    for piter = 1:npoints
        p.(par)                             = cs_grid(piter);
        [cs_sol{piter}, cs_stats{piter}]    = find_ss(opt, p, rates_old);
    end
    cs_res.cs_sol           = cs_sol;
    cs_res.cs_stats         = cs_stats;
    cs_res.cs_grid          = cs_grid;
    cs_res.par              = par;
end