function stats = calc_stats(opt, glob, p, s)
    TD              = 0;
    TB              = 0;
    TS              = 0;
    TDz             = zeros(p.Nz, 1);
    TBz             = zeros(p.Nz, 1);
    TSz             = zeros(p.Nz, 1);
    ldist           = cell(p.Nz, 1);
    ildist          = cell(p.Nz, 1);
    ldist_total     = zeros(p.Nb, 1);
    ildist_total    = zeros(1, p.Na);
    dist_vec        = reshape(s.dst, p.Nb * p.Na * p.Nz, 1);
    cons_vec        = reshape(s.cpol, p.Nb * p.Na * p.Nz, 1);
    h_vec           = reshape(s.hpol, p.Nb * p.Na * p.Nz, 1);
    dep_vec         = reshape(s.dpol, p.Nb * p.Na * p.Nz, 1);
    V_vec           = reshape(s.V, p.Nb * p.Na * p.Nz, 1);
    wt              = p.dtildea_vec .* p.dtildeb_vec .* dist_vec;
    wt_array        = reshape(wt, p.Nb, p.Na, p.Nz);
    z_vec           = reshape(p.zzz, p.Nb * p.Na * p.Nz, 1);
    H               = sum(wt .* z_vec .* h_vec);   % labor supply
    K               = H * p.K_H;        % get capital

  
    Rb              = p.r_plus .* (p.bbb > 0) + p.r_minus .* (p.bbb < 0);
    Ra              = p.r_F .* ones(p.Nb, p.Na, p.Nz);
    total_inc       = reshape(p.w * p.zzz .* s.hpol + Rb .* p.bbb + Ra .* p.aaa, p.Nb * p.Na * p.Nz, 1);
    
    %% Marginal distributions
    for nz = 1:p.Nz
        ldist{nz}       = trapz(p.a, s.dst(:, :, nz), 2);
        ildist{nz}      = trapz(p.b, s.dst(:, :, nz), 1);
        ldist_total     = ldist_total + ldist{nz};
        ildist_total    = ildist_total + ildist{nz};
    end

    
    %% Averages and totals
    TD              = sum(wt .* p.bfrombaz .* (p.bfrombaz > 0));
    TB              = sum(wt .* p.bfrombaz .* (p.bfrombaz < 0));
    TS              = sum(wt .* p.afrombaz);
    NW              = TS / (1 + p.mu_bank * (p.x_a - 1));                                              % Net worth = Illiquid - p.mu_bank * deposits 
    TD_bank         = (p.x_a - 1) * NW ;   
    
    
    a_mean          = TS;
    b_mean          = TD + TB;
    LIQ_ILLIQ       = b_mean / a_mean;
    
    cons_mean       = sum(wt .* cons_vec);
    total_inc_mean  = sum(wt .* total_inc);
    V_mean          = sum(wt .* V_vec);
    
    for nz = 1:p.Nz
        TDz(nz)     = sum(wt .* (p.zindfrombaz == nz) .* p.bfrombaz .* (p.bfrombaz > 0));
        TBz(nz)     = sum(wt .* (p.zindfrombaz == nz) .* p.bfrombaz .* (p.bfrombaz < 0));
        TSz(nz)     = sum(wt .* (p.zindfrombaz == nz) .* p.afrombaz);
    end
    
    %% Conditional group totals for the wealthy (w) and the poor (p)
    b_mean_w        = sum(wt .* p.bfrombaz .* (p.aindfrombaz > 1)) / sum(wt .* (p.aindfrombaz > 1));
    b_mean_p        = sum(wt .* p.bfrombaz .* (p.aindfrombaz == 1)) / sum(wt .* (p.aindfrombaz == 1));
    cons_mean_w     = sum(wt .* cons_vec .* (p.aindfrombaz > 1)) / sum(wt .* (p.aindfrombaz > 1));
    cons_mean_p     = sum(wt .* cons_vec .* (p.aindfrombaz == 1)) / sum(wt .* (p.aindfrombaz == 1));
    V_mean_w        = sum(wt .* V_vec .* (p.aindfrombaz > 1)) / sum(wt .* (p.aindfrombaz > 1));
    V_mean_p        = sum(wt .* V_vec .* (p.aindfrombaz == 1)) / sum(wt .* (p.aindfrombaz == 1));
    total_inc_mean_w    = sum(wt .* total_inc .* (p.aindfrombaz > 1)) / sum(wt .* (p.aindfrombaz > 1));
    total_inc_mean_p    = sum(wt .* total_inc .* (p.aindfrombaz == 1)) / sum(wt .* (p.aindfrombaz == 1));
    
    %% Inequality
    
    % Marginals
    b_marg                  = sum(wt_array, [2, 3]);
    a_marg                  = sum(wt_array, [1, 3]);
    z_marg                  = sum(wt_array, [1, 2]);
    
    [cons_vec_sorted, ord_c]  = sort(cons_vec);
    dcons_vec       = 0.5 * [cons_vec_sorted(1); 
                            cons_vec_sorted(3:end) - cons_vec_sorted(1:end - 2); 
                            cons_vec_sorted(end)];
    cons_marg_cum   = cumsum(wt(ord_c));
    cons_gini       = sum(cons_marg_cum .* (1 - cons_marg_cum) .* dcons_vec) / cons_mean;
    
    [a_vec_sorted, ord_a]     = sort(p.afrombaz);
    da_vec      	= 0.5 * [a_vec_sorted(1); 
                            a_vec_sorted(3:end) - a_vec_sorted(1:end - 2); 
                            a_vec_sorted(end)];
    a_marg_cum      = cumsum(wt(ord_a));
    a_gini          = sum(a_marg_cum .* (1 - a_marg_cum) .* da_vec) / TS;
    
    [b_vec_sorted, ord_b]     = sort(p.bfrombaz);
    db_vec      	= 0.5 * [b_vec_sorted(1); 
                            b_vec_sorted(3:end) - b_vec_sorted(1:end - 2); 
                            b_vec_sorted(end)];
    b_marg_cum      = cumsum(wt(ord_b));
    b_gini          = sum(b_marg_cum .* (1 - b_marg_cum) .* db_vec) / (TD + TB);
    
    a10             = prctilew(p.a, a_marg, 0.1);
    a90             = prctilew(p.a, a_marg, 0.9);
    b10             = prctilew(p.b, b_marg, 0.1);
    b90             = prctilew(p.b, b_marg, 0.9);
    c10             = prctilew(cons_vec, wt, 0.1);
    c90             = prctilew(cons_vec, wt, 0.9);
    b_var           = sum((p.bfrombaz - TD - TB) .^ 2 .* wt);
    a_var           = sum((p.afrombaz - TS) .^ 2 .* wt);
    cons_var        = sum((cons_vec - cons_mean) .^ 2 .* wt);
    
    
%     %% distributional IRFS
%     % need to align the grid and the distribution closer
%      [a_percentiles] = prctilew(p.a, a_marg, [0.1 0.25 0.5 0.75 0.99]);
%      [b_percentiles] = prctilew(p.b, b_marg, [0.1 0.25 0.5 0.75 0.9]);
%      
%      [~,ind_a]  = min( abs( p.a'-a_percentiles) );
%      % doing it in a dumb style
%     ginit               = zeros(p.Nb, p.Na, p.Nz);
%     ginit(:,ind_a(1),p.Nz) = 1;
%     ginit = reshape(ginit, p.Nb * p.Na * p.Nz, 1);
%     gtilde{1}               = reshape(ginit .* (p.dtildea_vec .* p.dtildeb_vec), p.Nb * p.Na, p.Nz);
%     gt{1}                   = reshape(ginit, p.Nb, p.Na, p.Nz);
% 
%      for t = 1:p.Nt-1
%             lmat                = eye(p.Nz, p.Nz) + p.dt(t) * p.la_mat_offdiag;
%         gtilde{t + 1}       = zeros(p.Nb * p.Na, p.Nz);
%     
%         for nz = 1:p.Nz
%             vec             = gtilde{t} * lmat(:, nz);
%             gtilde{t + 1}(:, nz)    = statst{t}.B{nz} \ vec;
%         end
%         gt{t + 1}           = reshape(gtilde{t + 1}, p.Nb * p.Na * p.Nz, 1);
%         gt{t + 1}           = reshape(gt{t + 1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
%     end
   

    %% Pack
    stats           = p;
    stats.dst       = s.dst;
    for nz = 1:p.Nz
        stats.ldist{nz}     = ldist{nz};
        stats.ildist{nz}    = ildist{nz};
    end
    stats.ldist_total   = ldist_total;
    stats.ildist_total  = ildist_total;
    stats.TD        = TD;
    stats.TB        = TB;
    stats.TS        = TS;
    stats.TDz       = TDz;
    stats.TBz       = TBz;
    stats.TSz       = TSz;
    stats.NW        = NW;
    stats.TD_bank   = TD_bank;
    stats.LIQ_ILLIQ = LIQ_ILLIQ;
    stats.cons_mean = cons_mean;
    stats.a_mean    = a_mean;
    stats.b_mean    = b_mean;
    stats.V_mean    = V_mean;
    stats.total_inc_mean    = total_inc_mean;
    stats.b_mean_w  = b_mean_w;
    stats.b_mean_p  = b_mean_p;
    stats.cons_mean_w   = cons_mean_w;
    stats.cons_mean_p   = cons_mean_p;
    stats.V_mean_w  = V_mean_w;
    stats.V_mean_p  = V_mean_p;
    stats.total_inc_mean_w  = total_inc_mean_w;
    stats.total_inc_mean_p  = total_inc_mean_p;
    stats.cons_var  = cons_var;
    stats.a_marg    = a_marg;
    stats.b_marg    = b_marg;
    stats.wt        = wt;
    stats.a_var     = a_var;
    stats.b_var     = b_var;
    stats.cons_gini = cons_gini;
    stats.a_gini    = a_gini;
    stats.b_gini    = b_gini;
    stats.c10       = c10;
    stats.c90       = c90;
    stats.a10       = a10;
    stats.a90       = a90;
    stats.b10       = b10;
    stats.b90       = b90;
    stats.K         = K;
    stats.H         = H;
    stats.r_plus    = p.r_plus;
    stats.r_minus   = p.r_minus;
    stats.r_F       = p.r_F;
    stats.inflow    = -(p.r_F) * NW; 
    stats.r_minus_div   = p.r_minus;
    stats.r_minus_gain  = 0;
    stats.cpol          = s.cpol;
    stats.hpol          = s.hpol;
    stats.dpol          = s.dpol;
end