function [guesses1, statst] = transition_iterate(opt, glob, p, guesses, init_state, final_ss)
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
    dt{p.Nt}        = final_ss.dpol;
    
    %% Guesses
    x_at            = guesses.x_at;
    Kt              = guesses.Kt;
    
    %% r_minus and w from K
    r_minust        = prod_K(opt, glob, p, Kt) - p.delta;
    wt              = prod_L(opt, glob, p, Kt);
    
    %% eta from leverage
    etat            = x_at * p.theta_bank;

    %% Use the differential equation for eta to get r_plus
    r_plust(p.Nt)   = final_ss.r_plus;
    for t = p.Nt-1:-1:1
        r_plust(t)  = (r_minust(t) * etat(t) ^ 2 / p.theta_bank - etat(t) * (p.rho_bank + p.f_bank) + p.f_bank + (etat(t + 1) - etat(t)) / p.dt(t)) / (etat(t) ^ 2 / p.theta_bank - etat(t));
    end
    
    %% Spread
    spreadt         = r_minust - r_plust;
    
    %% r_F
    r_Ft            = x_at .* (r_minust - r_plust) + r_plust;
    
    %% If GK then compute Nt and modify r_Ft
    if opt.GK
%         NWt(p.Nt)          = final_ss.NW;
%         for t = p.Nt-1:-1:1
%             NWt(t)      = (NWt(t+1) - p.MGK * p.dt(t)) / (1 + (x_at(t) .* (r_minust(t) - r_plust(t)) + r_plust(t) - p.f_bank) * p.dt(t));
%             r_Ft(t)     = p.f_bank - p.MGK / NWt(t);
%         end
%         NWt(1)      = final_ss.NW;
%         for t = 1:p.Nt-1
%             NWt(t + 1)  = NWt(t) * (1 + (x_at(t) .* (r_minust(t) - r_plust(t)) + r_plust(t) - p.f_bank) * p.dt(t)) + p.MGK * p.dt(t);
%             r_Ft(t)     = p.f_bank - p.MGK / NWt(t);
%         end
        % This implicit version is stable while the above are not
        NWt(1)      = init_state.NW;
        for t = 1:p.Nt-1
            NWt(t + 1)  = (NWt(t) + p.MGK * p.dt(t)) / (1 - (x_at(t+1) .* (r_minust(t+1) - r_plust(t+1)) + r_plust(t+1) - p.f_bank) * p.dt(t));
            r_Ft(t)     = p.f_bank - p.MGK / NWt(t);
        end
    end
    
    %% Solve HJB backward
    for t = p.Nt-1:-1:1
        p.r_plus                            = r_plust(t);
        p.r_minus                           = r_minust(t);
        p.r_F                               = r_Ft(t);
        if opt.GK
           p.r_F    = p.r_F * NWt(t);
        end
        p.w                                 = wt(t);
        [Vt{t}, dt{t}, ct{t}, ~, BU(t, :)]  = hjb_update(opt, glob, p, Vt{t + 1}, p.dt(t));
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
        pars.r_F            = r_Ft(t);
        if opt.GK
           pars.aaa         = p.aaa * NWt(t);
           pars.afrombaz    = p.afrombaz * NWt(t);
        end
        pars.w              = wt(t);
        
        statst{t}           = calc_stats(opt, glob, pars, struct('dst', gt{t}, 'cpol', ct{t}, 'dpol', dt{t}, 'V', Vt{t}));
        TSt(t)              = statst{t}.TS;
        TBt(t)              = statst{t}.TB;
        TDt(t)              = statst{t}.TD;
        statst{t}.r_plus    = r_plust(t);
        statst{t}.r_minus   = r_minust(t);
        statst{t}.spread    = spreadt(t);
        if opt.GK
            statst{t}.TS    = NWt(t);
        end
    end
    %% New guesses
    if opt.GK
        TSt                 = NWt;
    end
    
    Kt1                     = x_at .* TSt - TBt;
    guesses1.Kt             = opt.stepK * Kt1 + (1 - opt.stepK) * Kt;
    
    x_at1                   = TDt ./ TSt + 1;
    guesses1.x_at           = opt.stepK * x_at1 + (1 - opt.stepK) * x_at;
    
end