function guesses1 = transition_iterate2(opt, glob, p, guesses, init_state, final_ss)
    %% Preallocation
    BUt             = cell(p.Nt, p.Nz);
    Vt              = cell(p.Nt, 1);
    ct              = cell(p.Nt, 1);
    dt              = cell(p.Nt, 1);
    gt              = cell(p.Nt, 1);
    gtilde          = cell(p.Nt, 1);
    etat            = zeros(p.Nt, 1);
    TSt             = zeros(p.Nt, 1);
    TBt             = zeros(p.Nt, 1);
    TDt             = zeros(p.Nt, 1);
    
    %% Final period
    BUt(p.Nt, :)    = final_ss.BU;
    Vt{p.Nt}        = final_ss.V;
    ct{p.Nt}        = final_ss.cpol;
    dt{p.Nt}        = final_ss.dpol;
    
    %% Guesses
    r_plust         = guesses.r_plust;
    Kt              = guesses.Kt;
    
    %% r_minus and w from K
    r_minust        = prod_K(opt, glob, p, Kt) - p.delta;
    wt              = prod_L(opt, glob, p, Kt);
    spreadt         = r_minust - r_plust;
    % r_plust         = [r_plust; interp1(p.tgrid(1:end-1), r_plust, p.tgrid(end), 'linear', 'extrap')];
    % spreadt         = [spreadt; interp1(p.tgrid(1:end-1), spreadt, p.tgrid(end), 'linear', 'extrap')];
    
    %% Solve backwards the differential equation for eta with terminal condition from final_ss
%     [t, etat]       = ode45(@(t, x) min((p.rho_bank + p.f_bank - interp1(p.tgrid(1:end-1), r_plust, t, 'linear', 'extrap') ...
%                         - interp1(p.tgrid(1:end-1), spreadt, t, 'linear', 'extrap') .* x / p.theta_bank) .* x - p.f_bank, 0), ...
%                         flip(p.tgrid), final_ss.eta);
%     etat            = flip(etat);

    etat(p.Nt)      = final_ss.eta;
    for t = p.Nt-1:-1:1
        qA          = p.dt(t) * spreadt(t) / p.theta_bank;
        qB          = -(1 + p.dt(t) * (p.rho_bank + p.f_bank - r_plust(t)));
        qC          = p.dt(t) * p.f_bank + etat(t + 1);
        etat(t)     = (-qB - sqrt(qB ^ 2 - 4 * qA * qC)) / (2 * qA);
%         etat(t)     = etat(t + 1) - p.dt(t) * min((p.rho_bank + p.f_bank - r_plust(t + 1) - spreadt(t + 1) .* etat(t + 1) / p.theta_bank) .* etat(t + 1) - p.f_bank, 1e12);
    end
    
%     init_stats      = calc_stats(opt, glob, p, init_state);
%     etat(1)         = (init_stats.TS + init_stats.TD) / init_stats.TS * p.theta_bank;
%     for t = 2:p.Nt
% %         qA          = p.dt(t) * spreadt(t) / p.theta_bank;
% %         qB          = -(1 + p.dt(t) * (p.rho_bank + p.f_bank - r_plust(t)));
% %         qC          = p.dt(t) * p.f_bank + etat(t + 1);
% %         etat(t)     = (-qB - sqrt(qB ^ 2 - 4 * qA * qC)) / (2 * qA);
%          etat(t)    = etat(t - 1) + p.dt(t - 1) * min((p.rho_bank + p.f_bank - r_plust(t - 1) - spreadt(t - 1) .* etat(t - 1) / p.theta_bank) .* etat(t - 1) - p.f_bank, 1e12);
%     end
    
    %% Leverage
    x_at            = etat / p.theta_bank;
    
    %% r_F
    r_Ft            = x_at .* (r_minust - r_plust) + r_plust;
    
    %% Solve HJB backward
    for t = p.Nt-1:-1:1
        p.r_plus                            = r_plust(t);
        p.r_minus                           = r_minust(t);
        p.r_F                               = r_Ft(t);
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
            if t > -10
                vec         = gtilde{t} * lmat(:, nz);
            else
                vec         = gtilde0 * lmat(:, nz);
            end
            B                       = (1.0 - p.dt(t) * p.la_mat(nz, nz)) * speye(p.Nb * p.Na) - p.dt(t) * BU{t, nz}';

            gtilde{t + 1}(:, nz)    = B \ vec;
        end
        gt{t + 1}           = reshape(gtilde{t + 1}, p.Nb * p.Na * p.Nz, 1);
        gt{t + 1}           = reshape(gt{t + 1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
    end
    
    %% Get statistics
    for t = 1:p.Nt
        stats               = calc_stats(opt, glob, p, struct('dst', gt{t}, 'cpol', ct{t}, 'dpol', dt{t}, 'V', Vt{t}));
        TSt(t)              = stats.TS;
        TBt(t)              = stats.TB;
        TDt(t)              = stats.TD;
    end
    %% New guesses
    Kt1                     = x_at .* TSt - TBt;
    guesses1.Kt             = opt.stepK * Kt1 + (1 - opt.stepK) * Kt;
    
    TDe                     = TDt - (x_at .* TSt - TSt);
    % TDe                     = [TDe; (TDe(p.Nt) - TDe(p.Nt-1)) / p.dt(p.Nt-1) * p.dt(p.Nt) + TDe(p.Nt-1)];
    steprlevel              = 0.7 * (exp(-0.03 * p.tgrid(1:end)) - exp(-0.03 * p.tgrid(end))) + 0.05 + 0.1 * exp(-0.05 * (p.tgrid(end) - p.tgrid));
    % guesses1.r_plust        = [r_plust(1:end-1) - opt.stepr * 0.000001 * (TDe(2:end) - TDe(1:end-1)) ./ p.dt - opt.stepr * steprlevel .* TDe(1:end-1); final_ss.r_plus];
    
    % guesses1.r_plust        = r_plust - opt.stepr * steprlevel .* TDe;
    % guesses1.r_plust(1:end-1)   = guesses1.r_plust(1:end-1) + opt.stepr * 1 .* (TDe(2:end) - TDe(1:end-1)) ./ p.dt;
    
    guesses1.r_plust            = r_plust;
    guesses1.r_plust(1:end-1)   = r_plust(1:end-1) - opt.stepr * steprlevel(2:end) .* TDe(2:end) * 1.0;
    % guesses1.r_plust(2:end)     = guesses1.r_plust(2:end) + opt.stepr * steprlevel(2:end) .* TDe(1:end-1) * 0.1;
    guesses1.r_plust            = guesses1.r_plust + opt.stepr * steprlevel .* TDe * 0.05;
    
end