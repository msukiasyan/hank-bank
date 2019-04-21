function cs_res = comp_stat(opt, p, par, low, high, npoints)
    cs_grid                 = linspace(low, high, npoints);
    cs_sol                  = cell(npoints, 1);
    cs_stats                = cell(npoints, 1);
    rates_old               = [0.038 0.04];
    for piter = 1:npoints
        p.(par)             = cs_grid(piter);
        
        rates               = fmincon(@(x) eqs(opt, p, x), ...
                                rates_old, [1, -1], 0, [], [], [0.0 0.0], [0.1 0.1], [], ...
                                optimoptions('fmincon', 'Display', 'iter'));
        rates_old           = rates;

        p.r_F               = rates(2);
        p.r_minus           = rates(1);
        p                   = get_all_params(opt, p);
        
        cs_sol{piter}       = get_policies(opt, p);
        cs_stats{piter}     = calc_stats(opt, p, cs_sol{piter});
    end
    cs_res.cs_sol           = cs_sol;
    cs_res.cs_stats         = cs_stats;
end