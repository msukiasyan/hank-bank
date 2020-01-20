function [paths, statst]  = transition_Ashock_newton(opt, glob, p, sol, stats, sz, pers, guesses0)
    if nargin < 8
        guesses0.Kt         = ones(p.Nt, 1) * stats.K;
        guesses0.x_at       = repmat(stats.x_a, p.Nt, 1);
        guesses0.Nt         = ones(p.Nt, 1) * stats.N;
    end
    
    init_state          = sol;
    init_state.gvec     = reshape(init_state.dst, p.Nb * p.Na * p.Nz, 1);
    init_state.NW       = stats.NW;
    init_state.TB       = stats.TB;
    init_state.TD       = stats.TD;
    init_state.TD_bank  = stats.TD_bank;
    final_ss            = sol;
    final_ss.eta        = stats.eta;
    final_ss.r_plus     = stats.r_plus;
    final_ss.r_minus    = stats.r_minus;
    final_ss.r_F        = stats.r_F;
    final_ss.r_X        = stats.r_X;
    final_ss.spread     = stats.spread;
    final_ss.NW         = stats.NW;
    final_ss.K          = stats.K;
    final_ss.x_a        = stats.x_a;   
    final_ss.N          = stats.N;  
    if opt.GK
        p.NW                = 1;            % Make sure p.a is not distorted, because we want to 
                                            % lump changes in NW into r_F (more accurate)
        init_state.gvec     = init_state.gvec .* (stats.dtildea_vec .* stats.dtildeb_vec);
        p                   = setup(opt, glob, p);
        init_state.gvec     = init_state.gvec ./ (p.dtildea_vec .* p.dtildeb_vec);
    end
    
    p.Aprod                 = sz * exp(-pers * p.tgrid) + 1;
    
    [paths, statst]         = do_transition_newton(opt, glob, p, init_state, final_ss);
end