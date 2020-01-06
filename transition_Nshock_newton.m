function [paths, statst]  = transition_Nshock_newton(opt, glob, p, sol, stats, N0)
    guesses0.Kt         = ones(p.Nt, 1) * stats.K;
    guesses0.x_at       = repmat(stats.x_a, p.Nt, 1);
    
    init_state          = sol;
    for nb = 1:p.Nb
        for nz = 1:p.Nz
            init_state.dst(nb, :, nz)   = adjust_dist_prop(opt, glob, p, p.a, p.dtildea, init_state.dst(nb, :, nz), N0 / stats.NW);
        end
    end
    init_state.gvec     = reshape(init_state.dst, p.Nb * p.Na * p.Nz, 1);
    
    final_ss            = sol;
    final_ss.eta        = stats.eta;
    final_ss.r_plus     = stats.r_plus;
    final_ss.spread     = stats.spread;
    final_ss.NW         = stats.NW;
    final_ss.K          = stats.K;
    final_ss.x_a        = stats.x_a;
    
    if opt.GK
        p.NW                = 1;            % Make sure p.a is not distorted, because we want to 
                                            % lump changes in NW into r_F (more accurate)
        init_state.gvec     = init_state.gvec .* (stats.dtildea_vec .* stats.dtildeb_vec);
        p                   = setup(opt, glob, p);
        init_state.gvec     = init_state.gvec ./ (p.dtildea_vec .* p.dtildeb_vec);
    end
    
    [paths, statst]         = do_transition_newton(opt, glob, p, init_state, final_ss);
end