function [res, statst] = transition_residuals(opt, glob, p, guesses, init_state, final_ss)
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
    TD_bankt        = zeros(p.Nt, 1); 
    NWt             = zeros(p.Nt, 1);
    qt              = zeros(p.Nt, 1);
    Psit            = zeros(p.Nt, 1);
    iotat           = zeros(p.Nt, 1);
    r_minust        = zeros(p.Nt, 1);
    r_Xt            = zeros(p.Nt, 1);
  
    %% Final period
    BUt(p.Nt, :)    = final_ss.BU;
    Vt{p.Nt}        = final_ss.V;
    ct{p.Nt}        = final_ss.cpol;
    dt{p.Nt}        = final_ss.dpol;
    qt(p.Nt)        = 1;
    r_minust(p.Nt)  = final_ss.r_minus;
    
    %% Guesses
    x_at            = guesses.x_at;
    Kt              = guesses.Kt;
    
    %% Investment and capital price
    for t = 1:p.Nt-1
        Psit(t)     = (Kt(t + 1) - Kt(t)) / Kt(t) / p.dt(t) + p.delta;
        iotat(t)    = cap_prod_inv(opt, glob, p, Psit(t));
        qt(t)       = q_from_iota(opt, glob, p, iotat(t));
    end
    
    %% Adjust the initial distribution to reflect the change of q on impact
    N_change        = 1 - (1 - qt(1) / 1) * (init_state.NW + init_state.TD_bank - init_state.TB) / init_state.NW;
    for nb = 1:p.Nb
        for nz = 1:p.Nz
            init_state.dst(nb, :, nz)   = adjust_dist_prop(opt, glob, p, p.a, p.dtildea, init_state.dst(nb, :, nz), N_change);
        end
    end
    init_state.gvec     = reshape(init_state.dst, p.Nb * p.Na * p.Nz, 1);
    
    %% Return on capital
    mpk             = prod_K(opt, glob, p, Kt);
    for t = 1:p.Nt-1
        r_minust(t) = (mpk(t) - iotat(t) + (qt(t + 1) - qt(t)) / p.dt(t)) / qt(t) + Psit(t) - p.delta;
    end
    
    %% Wage from K
    wt              = prod_L(opt, glob, p, Kt);
    
    %% eta from leverage
    etat            = x_at * p.theta_bank;

    %% Use the differential equation for eta to get r_X
    r_Xt(p.Nt)      = final_ss.r_X;
    r_Ft(p.Nt)      = final_ss.r_F;
    r_plust(p.Nt)   = final_ss.r_plus;

    for t = p.Nt-1:-1:1
        r_Xt(t)     = (r_minust(t) * etat(t) ^ 2 / p.theta_bank - etat(t) * (p.rho_bank + p.f_bank) + p.f_bank + (etat(t + 1) - etat(t)) / p.dt(t)) / (etat(t) ^ 2 / p.theta_bank - etat(t));
        r_Ft(t)     = r_Xt(t) + (r_minust(t) - r_Xt(t)) * (etat(t) / p.theta_bank);
        r_plust(t)  = (r_Xt(t) - p.mu_bank * r_Ft(t)) / (1 - p.mu_bank); 
    end
    
    %% Spread
    spreadt         = r_minust - r_Xt;
    
    
    %% If GK then compute Nt and modify r_Ft
%     if opt.GK
% %         NWt(p.Nt)          = final_ss.NW;
% %         for t = p.Nt-1:-1:1
% %             NWt(t)      = (NWt(t+1) - p.MGK * p.dt(t)) / (1 + (x_at(t) .* (r_minust(t) - r_plust(t)) + r_plust(t) - p.f_bank) * p.dt(t));
% %             r_Ft(t)     = p.f_bank - p.MGK / NWt(t);
% %         end
% %         NWt(1)      = final_ss.NW;
% %         for t = 1:p.Nt-1
% %             NWt(t + 1)  = NWt(t) * (1 + (x_at(t) .* (r_minust(t) - r_plust(t)) + r_plust(t) - p.f_bank) * p.dt(t)) + p.MGK * p.dt(t);
% %             r_Ft(t)     = p.f_bank - p.MGK / NWt(t);
% %         end
%         % This implicit version is stable while the above are not
%         NWt(1)      = init_state.NW;
%         for t = 1:p.Nt-1
%             NWt(t + 1)  = (NWt(t) + p.MGK * p.dt(t)) / (1 - (x_at(t+1) .* (r_minust(t+1) - r_plust(t+1)) + r_plust(t+1) - p.f_bank) * p.dt(t));
%             r_Ft(t)     = p.f_bank - p.MGK / NWt(t);
%         end
%     end
    
    %% Calculate the transition and stats
    agg_paths.r_plust       = r_plust;
    agg_paths.r_minust      = r_minust;
    agg_paths.r_Ft          = r_Ft;
    agg_paths.r_Xt          = r_Xt;
    agg_paths.wt            = wt;
    agg_paths.x_at          = x_at;
    statst                  = transition_households(opt, glob, p, agg_paths, init_state, final_ss);
    Yt                      = p.Aprod .* Kt .^ p.alpha;
    
    for t = 1:p.Nt
        statst{t}.q         = qt(t);
        statst{t}.Y         = Yt(t);
    end
    
    %% residuals
    TSt                     = arrayfun(@(x) statst{x}.TS, (1:p.Nt)');
    TBt                     = arrayfun(@(x) statst{x}.TB, (1:p.Nt)');
    TDt                     = arrayfun(@(x) statst{x}.TD, (1:p.Nt)');
    
    NWt              =  TSt ./ (1 + p.mu_bank * (x_at - 1));                                              % Net worth = Illiquid - p.mu_bank * deposits 
    TD_bankt         =  TSt - NWt + TDt ;
    
     if opt.GK
        TSt                 = NWt;
    end
    
    Kt1                     = (x_at .* NWt - TBt) ./ qt;
    res.K                   = Kt1 ./ Kt - 1;
    
    x_at1                   = TD_bankt ./ NWt + 1;
    %x_at1                   = TD_bankt ./ TSt + 1;
    res.x_a                 = x_at1 ./ x_at - 1;
    
end