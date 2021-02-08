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
    NWt             = zeros(p.Nt, 1);
    qt              = zeros(p.Nt, 1);
    Psit            = zeros(p.Nt, 1);
    iotat           = zeros(p.Nt, 1);
    r_minust        = zeros(p.Nt, 1);
    
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
    piPt            = guesses.piPt;
    Ht              = guesses.Ht;
    
    %% Get inflation from 1) Taylor rule
    %% Get r_plus from 1) Taylor rule and 2) Fisher equation 
    r_plus_nomt = final_ss.r_plus + p.R_piP * piPt + p.MP; %% Taylor rule
    r_plust = r_plus_nomt - piPt; %% real rate = nominal rate - inflation
    
    %% Investment and capital price
    for t = 1:p.Nt-1
        Psit(t)     = (Kt(t + 1) - Kt(t)) / Kt(t) / p.dt(t) + p.delta;
        iotat(t)    = cap_prod_inv(opt, glob, p, Psit(t));
        qt(t)       = q_from_iota(opt, glob, p, iotat(t));
    end
    
    %% eta from leverage
    etat            = x_at * p.theta_bank;

    %% get output
    Yt                      = p.Aprod .* Kt .^ p.alpha .* Ht .^ (1 - p.alpha);

    %% use NKPC to get markups from inflation
    

    Ygrowtht      = [(Yt(2:end)-Yt(1:end-1))./(p.dt(1:end).*Yt(2:end));0];
    piPdott       = [(piPt(2:end)-piPt(1:end-1))./p.dt(1:end);0];
    invmarkupPt   = 1/p.markupP + p.sigP / p.zetaP * ((r_plust - Ygrowtht) .* piPt - piPdott);

    %invmarkupPt      = 1/p.markupP + p.sigP / p.zetaP * ((r_plust - Ygrowtht) .* piPt - piPdott);
    markupPt         = 1./invmarkupPt;
    %%  marginal product and factor prices
    mpk             = prod_K(opt, glob, p, Kt, Ht);
    wt              = prod_L(opt, glob, p, Kt, Ht) ./ markupPt;
    for t = 1:p.Nt-1
         r_minust(t) = (mpk(t) / markupPt(t) - iotat(t) + (qt(t + 1) - qt(t)) / p.dt(t)) / qt(t) + Psit(t) - p.delta;
        %r_minust(t) = (mpk(t) / markupPt(t) + (qt(t + 1) - qt(t)) / p.dt(t)) / qt(t) - p.delta;
    end
    
    r_plust1(p.Nt) = final_ss.r_plus;
    for t = p.Nt-1:-1:1
        r_plust1(t)  = (r_minust(t) * etat(t) ^ 2 / p.theta_bank - etat(t) * (p.rho_bank + p.f_bank) + p.f_bank + (etat(t + 1) - etat(t)) / p.dt(t)) / (etat(t) ^ 2 / p.theta_bank - etat(t));
    end
  
    %% get inflation implied by r_plust1
    piPt1 =  (r_plust1' - final_ss.r_plus - p.MP) / (p.R_piP-1);

    
    %% Spread
    spreadt         = r_minust - r_plust;
    
    %% r_F
    r_Ft            = x_at .* spreadt + r_plust;
    
    %% Adjust the initial distribution to reflect the change of q on impact
    N_change        = 1 - (1 - qt(1) / 1) * (init_state.NW + init_state.TD - init_state.TB) / init_state.NW;
    for nb = 1:p.Nb
        for nz = 1:p.Nz
            init_state.dst(nb, :, nz)   = adjust_dist_prop(opt, glob, p, p.a, p.dtildea, init_state.dst(nb, :, nz), N_change);
        end
    end
    init_state.gvec     = reshape(init_state.dst, p.Nb * p.Na * p.Nz, 1);
    


    
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
    
    %% Calculate the transition and stats
    agg_paths.r_plust       = r_plust;
    agg_paths.r_minust      = r_minust;
    agg_paths.r_Ft          = r_Ft;
    agg_paths.wt            = wt;
    agg_paths.x_at          = x_at;
    agg_paths.K_Ht          = Kt ./ Ht;
    
    statst                  = transition_households(opt, glob, p, agg_paths, init_state, final_ss);
    
    for t = 1:p.Nt
        statst{t}.r_plusnom = r_plus_nomt(t) ;
        statst{t}.piP       = piPt(t);
        statst{t}.q         = qt(t);
        statst{t}.Y         = Yt(t);
        statst{t}.I         = iotat(t) * Kt(t);
    end
    
    %% residuals
    TSt                     = arrayfun(@(x) statst{x}.TS, (1:p.Nt)');
    TBt                     = arrayfun(@(x) statst{x}.TB, (1:p.Nt)');
    TDt                     = arrayfun(@(x) statst{x}.TD, (1:p.Nt)');
    
    
    if opt.GK
        TSt                 = NWt;
    end
    
    Kt1                     = (x_at .* TSt - TBt) ./ qt;
    res.K                   = Kt1 ./ Kt - 1;
    
    x_at1                   = TDt ./ TSt + 1;
    res.x_a                 = x_at1 ./ x_at - 1;
    

    %res.piP                 = %exp(r_plust)./exp( r_plust1') - 1; 
    res.piP                = piPt1 -piPt; 
    Ht1                     = arrayfun(@(x) statst{x}.H, (1:p.Nt)');
    res.H                   = Ht1 ./ Ht - 1;
    
end