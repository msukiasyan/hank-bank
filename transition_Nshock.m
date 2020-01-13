function [paths, statst]  = transition_Nshock(opt, glob, p, sol, stats, N0)
    guesses0.Kt         = ones(p.Nt, 1) * stats.K;
    guesses0.x_at       = repmat(stats.x_a, p.Nt, 1);
    
    init_state          = sol;
    for nb = 1:p.Nb
        for nz = 1:p.Nz
            init_state.dst(nb, :, nz)   = adjust_dist_prop(opt, glob, p, p.a, p.dtildea, init_state.dst(nb, :, nz), N0 / stats.NW);
        end
    end
    init_state.gvec     = reshape(init_state.dst, p.Nb * p.Na * p.Nz, 1);
    init_state.NW       = N0;
    init_state.TB       = stats.TB;
    init_state.TD       = stats.TD;
    
    final_ss            = sol;
    final_ss.eta        = stats.eta;
    final_ss.r_plus     = stats.r_plus;
    final_ss.r_minus    = stats.r_minus;
    final_ss.spread     = stats.spread;
    final_ss.NW         = stats.NW;
    
    if opt.GK
        p.NW                = 1;            % Make sure p.a is not distorted, because we want to 
                                            % lump changes in NW into r_F (more accurate)
        init_state.gvec     = init_state.gvec .* (stats.dtildea_vec .* stats.dtildeb_vec);
        p                   = setup(opt, glob, p);
        init_state.gvec     = init_state.gvec ./ (p.dtildea_vec .* p.dtildeb_vec);
    end
    
    for it = 1:opt.maxittrans
        [guesses1, statst]  = transition_iterate(opt, glob, p, guesses0, init_state, final_ss);
        diffK               = max(abs(guesses1.Kt - guesses0.Kt));
        diffx_a             = max(abs(guesses1.x_at - guesses0.x_at));
        guesses0            = guesses1;
        disp(diffK);
        disp(diffx_a);
        if max(diffK, diffx_a) < opt.transtol
            break
        end
    end
    paths               = guesses0;
end