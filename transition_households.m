function statst = transition_households(opt, glob, p, agg_paths, init_state, final_ss)
    %% Preallocation
    statst          = cell(p.Nt, 1);
    BUt             = cell(p.Nt, p.Nz);
    Vt              = cell(p.Nt, 1);
    ct              = cell(p.Nt, 1);
    dt              = cell(p.Nt, 1);
    gt              = cell(p.Nt, 1);
    gtilde          = cell(p.Nt, 1);
    etat            = zeros(p.Nt, 1);
    r_plust         = zeros(p.Nt, 1);
    TSt             = zeros(p.Nt, 1);
    TBt             = zeros(p.Nt, 1);
    TDt             = zeros(p.Nt, 1);
    NWt             = zeros(p.Nt, 1);
    
    %% Final period
    BUt(p.Nt, :)    = final_ss.BU;
    Vt{p.Nt}        = final_ss.V;
    ct{p.Nt}        = final_ss.cpol;
    ht{p.Nt}        = final_ss.hpol;
    dt{p.Nt}        = final_ss.dpol;
    
    %% Aggregates
    r_plust         = agg_paths.r_plust;
    r_minust        = agg_paths.r_minust;
    x_at            = agg_paths.x_at;
    r_Ft            = agg_paths.r_Ft;
    wt              = agg_paths.wt;
    K_Ht            = agg_paths.K_Ht;
    %% Solve HJB backward
    for t = p.Nt-1:-1:1
        p.r_plus                            = r_plust(t);
        p.r_minus                           = r_minust(t);
        p.r_F                               = r_Ft(t);
        if opt.GK
           p.r_F    = p.r_F * NWt(t);
        end
        p.w                                 = wt(t);
        [Vt{t}, dt{t}, ct{t}, ht{t}, ~, BU(t, :)]  = hjb_update(opt, glob, p, Vt{t + 1}, p.dt(t));
    end

    %% Solve KFE forward
    gtilde{1}               = reshape(init_state.gvec .* (p.dtildea_vec .* p.dtildeb_vec), p.Nb * p.Na, p.Nz);
    gt{1}                   = reshape(init_state.gvec, p.Nb, p.Na, p.Nz);
    for t = 1:p.Nt-1
        lmat                = eye(p.Nz, p.Nz) + p.dt(t) * p.la_mat_offdiag;
        gtilde{t + 1}       = zeros(p.Nb * p.Na, p.Nz);
        for nz = 1:p.Nz
            vec             = gtilde{t} * lmat(:, nz);
            B               = (1.0 - p.dt(t) * p.la_mat(nz, nz)) * speye(p.Nb * p.Na) - p.dt(t) * BU{t, nz}';

            gtilde{t + 1}(:, nz)    = B \ vec;
        end
        gt{t + 1}           = reshape(gtilde{t + 1}, p.Nb * p.Na * p.Nz, 1);
        gt{t + 1}           = reshape(gt{t + 1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
    end
    
    %% Get statistics
    for t = 1:p.Nt
        pars                = p;
        pars.r_plus         = r_plust(t);
        pars.r_minus        = r_minust(t);
        pars.x_a            = x_at(t);
        pars.K_H            = K_Ht(t);
    
        if opt.GK
           pars.aaa         = p.aaa * NWt(t);
           pars.afrombaz    = p.afrombaz * NWt(t);
        end
        pars.w              = wt(t);
        
        statst{t}           = calc_stats(opt, glob, pars, struct('dst', gt{t}, 'cpol', ct{t}, 'hpol', ht{t}, 'dpol', dt{t}, 'V', Vt{t}));
        TSt(t)              = statst{t}.TS;
        TBt(t)              = statst{t}.TB;
        TDt(t)              = statst{t}.TD;
        TDt(t)              = statst{t}.TD;
        statst{t}.r_plus    = r_plust(t);
        statst{t}.r_minus   = r_minust(t);
        statst{t}.r_F       = r_Ft(t);
        statst{t}.w         = wt(t);
        statst{t}.spread    = r_minust(t) - r_plust(t);
        if opt.GK
            statst{t}.TS    = NWt(t);
        end
    end
    
end