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
    
    %% Calculate the transition and stats
    agg_paths.r_plust       = r_plust;
    agg_paths.r_minust      = r_minust;
    agg_paths.r_Ft          = r_Ft;
    agg_paths.wt            = wt;
    
    statst                  = transition_households(opt, glob, p, agg_paths, init_state, final_ss);
    
    %% residuals
    TSt                     = arrayfun(@(x) statst{x}.TS, (1:p.Nt)');
    TBt                     = arrayfun(@(x) statst{x}.TB, (1:p.Nt)');
    TDt                     = arrayfun(@(x) statst{x}.TD, (1:p.Nt)');
    
    if opt.GK
        TSt                 = NWt;
    end
    
    Kt1                     = x_at .* TSt - TBt;
    res.K                   = Kt1 ./ Kt - 1;
    
    x_at1                   = TDt ./ TSt + 1;
    res.x_a                 = x_at1 ./ x_at - 1;
    
end