function paths  = transition_Nshock2(opt, p, sol, stats, N0)
    guesses0.Kt         = ones(p.Nt, 1) * stats.K;
    guesses0.r_plust    = repmat(stats.r_plus, p.Nt, 1);
    
    init_state          = sol;
    for nb = 1:p.Nb
        for nz = 1:p.Nz
            init_state.dst(nb, :, nz)   = adjust_dist_prop(opt, p, p.a, p.dtildea, init_state.dst(nb, :, nz), N0 / stats.NW);
        end
    end
    init_state.gvec     = reshape(init_state.dst, p.Nb * p.Na * p.Nz, 1);
    
    final_ss            = sol;
    final_ss.eta        = stats.eta;
    final_ss.r_plus     = stats.r_plus;
    final_ss.spread     = stats.spread;
    for it = 1:opt.maxittrans
        guesses1        = transition_iterate2(opt, p, guesses0, init_state, final_ss);
        diffK           = max(abs(guesses1.Kt - guesses0.Kt));
        diffr_plus      = max(abs(guesses1.r_plust - guesses0.r_plust));
        guesses0        = guesses1;
        disp(diffK);
        disp(diffr_plus);
        if max(diffK, diffr_plus) < opt.transtol
            break
        end
    end
    paths               = guesses0;
end