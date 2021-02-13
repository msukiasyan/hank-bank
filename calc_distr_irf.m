function [statst] = calc_distr_irf(opt, glob, p,stats, statst, statst_benchmark, weighted)


    
     %% distributional IRFS
     %% first do it by illiquid percentiles
     a_vec = p.afrombaz; 
     b_vec = p.bfrombaz; 
    % need to align the grid and the distribution closer
     [a_percentiles_pre] = prctilew(p.a, statst_benchmark{1}.a_marg, [0.0001 0.1 0.25 0.5 0.75 0.9 0.95]);
     [b_percentiles]     = prctilew(p.b, statst_benchmark{1}.b_marg, [0.0001 0.1 0.25 0.5 0.75 0.9 0.95]);
        
     % adjust percentiles for the initial jump
     [a_percentiles] =  (statst{1}.TS /statst_benchmark{1}.TS) * a_percentiles_pre;
     % [a_percentiles] =  prctilew(p.a, statst{1}.a_marg, [0.0001 0.1 0.25 0.5 0.75 0.9 0.95]);
        
     % correct for the 0th 
     a_percentiles_pre(1) =0;
     a_percentiles(1)     =0;
     
     % doing it in a dumb way
     a_rep          = repmat(p.a,7,1);
     a_percentiles_rep = repmat(a_percentiles',1,p.Na);
     a_repNaN       = a_rep;
     a_repNaN(a_rep<a_percentiles_rep) = NaN;
     
     a_percentiles_pre_rep = repmat(a_percentiles_pre',1,p.Na);
     a_pre_repNaN       = a_rep;
     a_pre_repNaN(a_rep<a_percentiles_pre_rep) = NaN;
     
     
     [~,i]          = min(a_repNaN,[],2);
     a_ind_right    = i;
     a_ind_left     = max(1,i-1);
     a_weight_left  = 1 - (a_percentiles - a_rep(1,a_ind_left))./(a_rep(1,a_ind_right) - a_rep(1,a_ind_left));
     a_weight_right = 1 - a_weight_left;
     a_weight_left(i==1) = 1;
     a_weight_right(i==1) = 0;
     
     [~,j]          = min(a_pre_repNaN,[],2);
     a_pre_ind_right    = j;
     a_pre_ind_left     = max(1,j-1);
     a_pre_weight_left  = 1- (a_percentiles_pre - a_rep(1,a_pre_ind_left))./(a_rep(1,a_pre_ind_right) - a_rep(1,a_pre_ind_left));
     a_pre_weight_right = 1 - a_pre_weight_left; 
     a_pre_weight_left(j==1) = 1;
     a_pre_weight_right(j==1) = 0;
     
     if weighted == 0
     a_weight_left(a_weight_left>=0.5) = 1;
     a_weight_right(a_weight_right<=0.5) = 0;
     a_pre_weight_left(a_pre_weight_left>=0.5) = 1;
     a_pre_weight_right(a_pre_weight_right<=0.5) = 0;
     end
     


     for k = 1:7

    cons_pre_vec            = reshape(statst_benchmark{1}.cpol, p.Nb * p.Na * p.Nz, 1);
    h_pre_vec               = reshape(statst_benchmark{1}.hpol, p.Nb * p.Na * p.Nz, 1);
    dep_pre_vec             = reshape(statst_benchmark{1}.dpol, p.Nb * p.Na * p.Nz, 1);
    s_pre_vec               = reshape(statst_benchmark{1}.spol, p.Nb * p.Na * p.Nz, 1);
    m_pre_vec               = reshape(statst_benchmark{1}.mpol, p.Nb * p.Na * p.Nz, 1);
    labor_inc_pre_vec       = reshape(statst_benchmark{1}.labor_inc, p.Nb * p.Na * p.Nz, 1);
    liq_inc_pre_vec         = reshape(statst_benchmark{1}.liq_inc, p.Nb * p.Na * p.Nz, 1);
    illiq_inc_pre_vec       = reshape(statst_benchmark{1}.illiq_inc, p.Nb * p.Na * p.Nz, 1);
    disposable_inc_pre_vec        = reshape(statst_benchmark{1}.disposable_inc, p.Nb * p.Na * p.Nz, 1);


    
    cons_vec                = reshape(statst{1}.cpol, p.Nb * p.Na * p.Nz, 1);
    h_vec               	= reshape(statst{1}.hpol, p.Nb * p.Na * p.Nz, 1);
    dep_vec                 = reshape(statst{1}.dpol, p.Nb * p.Na * p.Nz, 1);
    s_vec                   = reshape(statst{1}.spol, p.Nb * p.Na * p.Nz, 1);
    m_vec               	= reshape(statst{1}.mpol, p.Nb * p.Na * p.Nz, 1);
    labor_inc_vec           = reshape(statst{1}.labor_inc, p.Nb * p.Na * p.Nz, 1);
    liq_inc_vec             = reshape(statst{1}.liq_inc, p.Nb * p.Na * p.Nz, 1);
    illiq_inc_vec           = reshape(statst{1}.illiq_inc, p.Nb * p.Na * p.Nz, 1);    
    disposable_inc_vec      = reshape(statst{1}.disposable_inc, p.Nb * p.Na * p.Nz, 1);
 
    ginit_pre_left     =  zeros(p.Nb, p.Na, p.Nz);
     ginit_pre_right    =  zeros(p.Nb, p.Na, p.Nz);
        
     ginit_left         =  zeros(p.Nb, p.Na, p.Nz);
     ginit_right        =  zeros(p.Nb, p.Na, p.Nz);
     
     ginit_left(:,a_ind_left(k),:)   =   statst{1}.dst(:,a_ind_left(k),:);
     ginit_right(:,a_ind_right(k),:) =   statst{1}.dst(:,a_ind_right(k),:);
     
     %ginit_pre_left(:,a_pre_ind_left(k),:)   =   statst_benchmark{1}.dst(:,a_pre_ind_left(k),:);
     %ginit_pre_right(:,a_pre_ind_right(k),:) =   statst_benchmark{1}.dst(:,a_pre_ind_right(k),:);
     
         
     ginit_pre_left(:,a_pre_ind_left(k),:)   =   stats.dst(:,a_pre_ind_left(k),:);
     ginit_pre_right(:,a_pre_ind_right(k),:) =   stats.dst(:,a_pre_ind_right(k),:);
     
     
     ginit                           =   a_weight_left(k) .* ginit_left + a_weight_right(k) .* ginit_right;
     ginit_pre                       =   a_pre_weight_left(k) .* ginit_pre_left + a_pre_weight_right(k) .* ginit_pre_right;
    
     % do the impact adjustment for the people with 0 illiquid assets
     if k == 1
     ginit                           = ginit_pre;
     end
     
     ginit                   = reshape(ginit, p.Nb * p.Na * p.Nz, 1); 
     gtilde{1}               = reshape(ginit .* (p.dtildea_vec .* p.dtildeb_vec), p.Nb * p.Na, p.Nz);
     gt{1}                   = reshape(gtilde{1}, p.Nb * p.Na * p.Nz, 1);
     gt{1}                   = reshape(gt{1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
     dist_vec                = reshape(gt{1}, p.Nb * p.Na * p.Nz, 1);
     wt{1}                   = p.dtildea_vec .* p.dtildeb_vec .*  dist_vec  ;
     wt{1}                   = wt{1} / sum(wt{1}); 
    
     ginit_pre               = reshape(ginit_pre, p.Nb * p.Na * p.Nz, 1); 
     gtilde_pre{1}           = reshape(ginit_pre .* (p.dtildea_vec .* p.dtildeb_vec), p.Nb * p.Na, p.Nz);
     gt_pre{1}               = reshape(gtilde_pre{1}, p.Nb * p.Na * p.Nz, 1);
     gt_pre{1}               = reshape(gt_pre{1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
     dist_pre_vec            = reshape(gt_pre{1}, p.Nb * p.Na * p.Nz, 1);
     wt_pre{1}               = p.dtildea_vec .* p.dtildeb_vec .*  dist_pre_vec  ;
     wt_pre{1}               = wt_pre{1} / sum(wt_pre{1}); 
    


    
    statst{1}.illiq.cons_irf_percentile(k)        = (sum(wt{1} .* cons_vec) - sum(wt_pre{1} .* cons_pre_vec)) ./ sum(wt_pre{1} .* cons_pre_vec) * 100;
    statst{1}.illiq.h_irf_percentile(k)           = (sum(wt{1} .* h_vec) - sum(wt_pre{1} .* h_pre_vec)) / sum(wt_pre{1} .* h_pre_vec) * 100;
    statst{1}.illiq.dep_irf_percentile(k)         = (sum(wt{1} .* dep_vec) - sum(wt_pre{1} .* dep_pre_vec)); 
    statst{1}.illiq.s_irf_percentile(k)           = (sum(wt{1} .* s_vec) - sum(wt_pre{1} .* s_pre_vec));
    statst{1}.illiq.m_irf_percentile(k)           = (sum(wt{1} .* m_vec) - sum(wt_pre{1} .* m_pre_vec)) ;
    statst{1}.illiq.labor_inc_irf_percentile(k)   = (sum(wt{1} .* labor_inc_vec) - sum(wt_pre{1} .* labor_inc_pre_vec)) / sum(wt_pre{1} .* labor_inc_pre_vec) * 100;    
    statst{1}.illiq.liq_inc_irf_percentile(k)     = (sum(wt{1} .* liq_inc_vec) - sum(wt_pre{1} .* liq_inc_pre_vec)) ;
    statst{1}.illiq.illiq_inc_irf_percentile(k)   = (sum(wt{1} .* illiq_inc_vec) - sum(wt_pre{1} .* illiq_inc_pre_vec));    
    statst{1}.illiq.a_irf_percentile(k)           = (sum(wt{1} .* a_vec) - sum(wt_pre{1} .* a_vec));
    statst{1}.illiq.b_irf_percentile(k)           = (sum(wt{1} .* b_vec) - sum(wt_pre{1} .* b_vec));    
    statst{1}.illiq.disposable_inc_irf_percentile(k)   = (sum(wt{1} .* disposable_inc_vec) - sum(wt_pre{1} .* disposable_inc_pre_vec))  / sum(wt_pre{1} .*  disposable_inc_pre_vec) * 100;;  

    
    
    for t = 1:p.Nt-1
        lmat                    = eye(p.Nz, p.Nz) + p.dt(t) * p.la_mat_offdiag;
        gtilde{t + 1}           = zeros(p.Nb * p.Na, p.Nz);
        gtilde_pre{t + 1}       = zeros(p.Nb * p.Na, p.Nz);

        for nz = 1:p.Nz
            vec             = gtilde{t} * lmat(:, nz);
            vec_pre         = gtilde_pre{t} * lmat(:, nz);
            gtilde{t + 1}(:, nz)        = statst{t}.B{nz} \ vec;
            gtilde_pre{t + 1}(:, nz)    = statst_benchmark{t}.B{nz} \ vec_pre;
        
        end
        
        gt{t + 1}           = reshape(gtilde{t + 1}, p.Nb * p.Na * p.Nz, 1);
        gt{t + 1}           = reshape(gt{t + 1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);  
        dist_vec            = reshape(gt{t+1}, p.Nb * p.Na * p.Nz, 1);
        wt{t + 1}           = p.dtildea_vec .* p.dtildeb_vec .* dist_vec;
        wt{t + 1}           = wt{t+1} / sum(wt{t+1}); 
        
        
        gt_pre{t + 1}           = reshape(gtilde_pre{t + 1}, p.Nb * p.Na * p.Nz, 1);
        gt_pre{t + 1}           = reshape(gt_pre{t + 1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
        dist_pre_vec            = reshape(gt_pre{t+1}, p.Nb * p.Na * p.Nz, 1);
        wt_pre{t + 1}           = p.dtildea_vec .* p.dtildeb_vec .* dist_pre_vec;
        wt_pre{t + 1}           = wt_pre{t+1} / sum(wt_pre{t+1});  
        
        cons_vec            = reshape(statst{t+1}.cpol, p.Nb * p.Na * p.Nz, 1);
        h_vec               = reshape(statst{t+1}.hpol, p.Nb * p.Na * p.Nz, 1);
        dep_vec             = reshape(statst{t+1}.dpol, p.Nb * p.Na * p.Nz, 1);
        s_vec               = reshape(statst{t+1}.spol, p.Nb * p.Na * p.Nz, 1);
        m_vec               = reshape(statst{t+1}.mpol, p.Nb * p.Na * p.Nz, 1);
        labor_inc_vec       = reshape(statst{t+1}.labor_inc, p.Nb * p.Na * p.Nz, 1);
        liq_inc_vec         = reshape(statst{t+1}.liq_inc, p.Nb * p.Na * p.Nz, 1);
        illiq_inc_vec       = reshape(statst{t+1}.illiq_inc, p.Nb * p.Na * p.Nz, 1);
        disposable_inc_vec  = reshape(statst{t+1}.disposable_inc, p.Nb * p.Na * p.Nz, 1);

        
        cons_pre_vec            = reshape(statst_benchmark{t+1}.cpol, p.Nb * p.Na * p.Nz, 1);
        h_pre_vec               = reshape(statst_benchmark{t+1}.hpol, p.Nb * p.Na * p.Nz, 1);
        dep_pre_vec             = reshape(statst_benchmark{t+1}.dpol, p.Nb * p.Na * p.Nz, 1);
        s_pre_vec               = reshape(statst_benchmark{t+1}.spol, p.Nb * p.Na * p.Nz, 1);
        m_pre_vec               = reshape(statst_benchmark{t+1}.mpol, p.Nb * p.Na * p.Nz, 1);
        labor_inc_pre_vec       = reshape(statst_benchmark{t+1}.labor_inc, p.Nb * p.Na * p.Nz, 1);
        liq_inc_pre_vec         = reshape(statst_benchmark{t+1}.liq_inc, p.Nb * p.Na * p.Nz, 1);
        illiq_inc_pre_vec       = reshape(statst_benchmark{t+1}.illiq_inc, p.Nb * p.Na * p.Nz, 1);
        disposable_inc_pre_vec  = reshape(statst_benchmark{t+1}.disposable_inc, p.Nb * p.Na * p.Nz, 1);

        
        statst{t+1}.illiq.cons_irf_percentile(k)          = (sum(wt{t + 1} .* cons_vec) - sum(wt_pre{t + 1} .* cons_pre_vec)) ./ sum(wt_pre{t + 1} .* cons_pre_vec) * 100;
        statst{t+1}.illiq.h_irf_percentile(k)             = (sum(wt{t + 1} .* h_vec) - sum(wt_pre{t + 1} .* h_pre_vec)) ./ sum(wt_pre{t + 1} .* h_pre_vec) * 100;
        statst{t+1}.illiq.dep_irf_percentile(k)           = (sum(wt{t + 1} .* dep_vec) - sum(wt_pre{t + 1} .* dep_pre_vec));
        statst{t+1}.illiq.s_irf_percentile(k)             = (sum(wt{t + 1} .* s_vec) - sum(wt_pre{t + 1} .* s_pre_vec)) ;
        statst{t+1}.illiq.m_irf_percentile(k)             = (sum(wt{t + 1} .* m_vec) - sum(wt_pre{t + 1} .* m_pre_vec)) ;
        statst{t+1}.illiq.labor_inc_irf_percentile(k)     = (sum(wt{t + 1} .* labor_inc_vec) - sum(wt_pre{t + 1} .* labor_inc_pre_vec)) ./ sum(wt_pre{t + 1} .* labor_inc_pre_vec) * 100;
        statst{t+1}.illiq.liq_inc_irf_percentile(k)       = (sum(wt{t + 1} .* liq_inc_vec) - sum(wt_pre{t + 1} .* liq_inc_pre_vec)) ;
        statst{t+1}.illiq.illiq_inc_irf_percentile(k)     = (sum(wt{t + 1} .* illiq_inc_vec) - sum(wt_pre{t + 1} .* illiq_inc_pre_vec));
        statst{t+1}.illiq.a_irf_percentile(k)             = (sum(wt{t + 1} .* a_vec) - sum(wt_pre{t + 1} .* a_vec)) ;
        statst{t+1}.illiq.b_irf_percentile(k)             = (sum(wt{t + 1} .* b_vec) - sum(wt_pre{t + 1} .* b_vec));
        statst{t+1}.illiq.disposable_inc_irf_percentile(k)     = (sum(wt{t + 1} .* disposable_inc_vec) - sum(wt_pre{t + 1} .* disposable_inc_pre_vec)) ./ sum(wt_pre{t + 1} .* disposable_inc_pre_vec) * 100;

    end
    
    end
   
   
    %% now repeat the exercise but with repect to liquid percentiles
     a_vec = p.afrombaz; 
     b_vec = p.bfrombaz; 
    % need to align the grid and the distribution closer
     [b_percentiles]     = prctilew(p.b, statst_benchmark{1}.b_marg, [0.0001 0.1 0.25 0.5 0.75 0.9 0.95]);
      
        
     % correct for the 0th 
     b_percentiles(1)     =0;
     b_percentiles_pre = b_percentiles;
     % doing it in a dumb way
     b_rep          = repmat(p.b',7,1);
     b_percentiles_rep = repmat(b_percentiles',1,p.Nb);
     b_repNaN       = b_rep;
     b_repNaN(b_rep<b_percentiles_rep) = NaN;
     
     b_percentiles_pre_rep = repmat(b_percentiles_pre',1,p.Nb);
     b_pre_repNaN       = b_rep;
     b_pre_repNaN(b_rep<b_percentiles_pre_rep) = NaN;
     
     
     [~,i]          = min(b_repNaN,[],2);
     b_ind_right    = i;
     b_ind_left     = max(1,i-1);
     b_weight_left  = 1 - (b_percentiles - b_rep(1,b_ind_left))./(b_rep(1,b_ind_right) - b_rep(1,b_ind_left));
     b_weight_right = 1 - b_weight_left;
     b_weight_left(i==1) = 1;
     b_weight_right(i==1) = 0;
     
     [~,j]          = min(b_pre_repNaN,[],2);
     b_pre_ind_right    = j;
     b_pre_ind_left     = max(1,j-1);
     b_pre_weight_left  = 1- (b_percentiles_pre - b_rep(1,b_pre_ind_left))./(b_rep(1,b_pre_ind_right) - b_rep(1,b_pre_ind_left));
     b_pre_weight_right = 1 - b_pre_weight_left; 
     b_pre_weight_left(j==1) = 1;
     b_pre_weight_right(j==1) = 0;
     
     if weighted == 0
     b_weight_left(b_weight_left>=0.5) = 1;
     b_weight_right(b_weight_right<=0.5) = 0;
     b_pre_weight_left(b_pre_weight_left>=0.5) = 1;
     b_pre_weight_right(b_pre_weight_right<=0.5) = 0;
     end
     


     for k = 1:7

    cons_pre_vec            = reshape(statst_benchmark{1}.cpol, p.Nb * p.Na * p.Nz, 1);
    h_pre_vec               = reshape(statst_benchmark{1}.hpol, p.Nb * p.Na * p.Nz, 1);
    dep_pre_vec             = reshape(statst_benchmark{1}.dpol, p.Nb * p.Na * p.Nz, 1);
    s_pre_vec               = reshape(statst_benchmark{1}.spol, p.Nb * p.Na * p.Nz, 1);
    m_pre_vec               = reshape(statst_benchmark{1}.mpol, p.Nb * p.Na * p.Nz, 1);
    labor_inc_pre_vec       = reshape(statst_benchmark{1}.labor_inc, p.Nb * p.Na * p.Nz, 1);
    liq_inc_pre_vec         = reshape(statst_benchmark{1}.liq_inc, p.Nb * p.Na * p.Nz, 1);
    illiq_inc_pre_vec       = reshape(statst_benchmark{1}.illiq_inc, p.Nb * p.Na * p.Nz, 1);
    disposable_inc_pre_vec  = reshape(statst_benchmark{1}.disposable_inc, p.Nb * p.Na * p.Nz, 1);


    
    cons_vec                = reshape(statst{1}.cpol, p.Nb * p.Na * p.Nz, 1);
    h_vec               	= reshape(statst{1}.hpol, p.Nb * p.Na * p.Nz, 1);
    dep_vec                 = reshape(statst{1}.dpol, p.Nb * p.Na * p.Nz, 1);
    s_vec                   = reshape(statst{1}.spol, p.Nb * p.Na * p.Nz, 1);
    m_vec               	= reshape(statst{1}.mpol, p.Nb * p.Na * p.Nz, 1);
    labor_inc_vec           = reshape(statst{1}.labor_inc, p.Nb * p.Na * p.Nz, 1);
    liq_inc_vec             = reshape(statst{1}.liq_inc, p.Nb * p.Na * p.Nz, 1);
    illiq_inc_vec           = reshape(statst{1}.illiq_inc, p.Nb * p.Na * p.Nz, 1);    
    disposable_inc_vec      = reshape(statst{1}.disposable_inc, p.Nb * p.Na * p.Nz, 1);    

     ginit_pre_left     =  zeros(p.Nb, p.Na, p.Nz);
     ginit_pre_right    =  zeros(p.Nb, p.Na, p.Nz);
        
     ginit_left         =  zeros(p.Nb, p.Na, p.Nz);
     ginit_right        =  zeros(p.Nb, p.Na, p.Nz);
     
     ginit_left(b_ind_left(k),:,:)   =   statst{1}.dst(b_ind_left(k),:,:);
     ginit_right(b_ind_right(k),:,:) =   statst{1}.dst(b_ind_right(k),:,:);
     
     %ginit_pre_left(:,a_pre_ind_left(k),:)   =   statst_benchmark{1}.dst(:,a_pre_ind_left(k),:);
     %ginit_pre_right(:,a_pre_ind_right(k),:) =   statst_benchmark{1}.dst(:,a_pre_ind_right(k),:);
     
         
     ginit_pre_left(b_pre_ind_left(k),:,:)   =   stats.dst(b_pre_ind_left(k),:,:);
     ginit_pre_right(b_pre_ind_right(k),:,:) =   stats.dst(b_pre_ind_right(k),:,:);
     
     
     ginit                           =   b_weight_left(k) .* ginit_left + b_weight_right(k) .* ginit_right;
     ginit_pre                       =   b_pre_weight_left(k) .* ginit_pre_left + b_pre_weight_right(k) .* ginit_pre_right;
    
     % do the impact adjustment for the people with 0 illiquid assets
     if k == 1
     ginit                           = ginit_pre;
     end
     
     ginit                   = reshape(ginit, p.Nb * p.Na * p.Nz, 1); 
     gtilde{1}               = reshape(ginit .* (p.dtildea_vec .* p.dtildeb_vec), p.Nb * p.Na, p.Nz);
     gt{1}                   = reshape(gtilde{1}, p.Nb * p.Na * p.Nz, 1);
     gt{1}                   = reshape(gt{1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
     dist_vec                = reshape(gt{1}, p.Nb * p.Na * p.Nz, 1);
     wt{1}                   = p.dtildea_vec .* p.dtildeb_vec .*  dist_vec  ;
     wt{1}                   = wt{1} / sum(wt{1}); 
    
     ginit_pre               = reshape(ginit_pre, p.Nb * p.Na * p.Nz, 1); 
     gtilde_pre{1}           = reshape(ginit_pre .* (p.dtildea_vec .* p.dtildeb_vec), p.Nb * p.Na, p.Nz);
     gt_pre{1}               = reshape(gtilde_pre{1}, p.Nb * p.Na * p.Nz, 1);
     gt_pre{1}               = reshape(gt_pre{1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
     dist_pre_vec            = reshape(gt_pre{1}, p.Nb * p.Na * p.Nz, 1);
     wt_pre{1}               = p.dtildea_vec .* p.dtildeb_vec .*  dist_pre_vec  ;
     wt_pre{1}               = wt_pre{1} / sum(wt_pre{1}); 
    


    
    statst{1}.liq.cons_irf_percentile(k)        = (sum(wt{1} .* cons_vec) - sum(wt_pre{1} .* cons_pre_vec)) ./ sum(wt_pre{1} .* cons_pre_vec) * 100;
    statst{1}.liq.h_irf_percentile(k)           = (sum(wt{1} .* h_vec) - sum(wt_pre{1} .* h_pre_vec)) / sum(wt_pre{1} .* h_pre_vec) * 100;
    statst{1}.liq.dep_irf_percentile(k)         = (sum(wt{1} .* dep_vec) - sum(wt_pre{1} .* dep_pre_vec)); 
    statst{1}.liq.s_irf_percentile(k)           = (sum(wt{1} .* s_vec) - sum(wt_pre{1} .* s_pre_vec));
    statst{1}.liq.m_irf_percentile(k)           = (sum(wt{1} .* m_vec) - sum(wt_pre{1} .* m_pre_vec)) ;
    statst{1}.liq.labor_inc_irf_percentile(k)   = (sum(wt{1} .* labor_inc_vec) - sum(wt_pre{1} .* labor_inc_pre_vec)) / sum(wt_pre{1} .* labor_inc_pre_vec) * 100;    
    statst{1}.liq.liq_inc_irf_percentile(k)     = (sum(wt{1} .* liq_inc_vec) - sum(wt_pre{1} .* liq_inc_pre_vec)) ;
    statst{1}.liq.illiq_inc_irf_percentile(k)   = (sum(wt{1} .* illiq_inc_vec) - sum(wt_pre{1} .* illiq_inc_pre_vec));    
    statst{1}.liq.a_irf_percentile(k)           = (sum(wt{1} .* a_vec) - sum(wt_pre{1} .* a_vec));
    statst{1}.liq.b_irf_percentile(k)           = (sum(wt{1} .* b_vec) - sum(wt_pre{1} .* b_vec));    
    statst{1}.liq.disposable_inc_irf_percentile(k)   = (sum(wt{1} .* disposable_inc_vec) - sum(wt_pre{1} .* disposable_inc_pre_vec)) / sum(wt_pre{1} .* disposable_inc_pre_vec) * 100;      

    
    
    for t = 1:p.Nt-1
        lmat                    = eye(p.Nz, p.Nz) + p.dt(t) * p.la_mat_offdiag;
        gtilde{t + 1}           = zeros(p.Nb * p.Na, p.Nz);
        gtilde_pre{t + 1}       = zeros(p.Nb * p.Na, p.Nz);

        for nz = 1:p.Nz
            vec             = gtilde{t} * lmat(:, nz);
            vec_pre         = gtilde_pre{t} * lmat(:, nz);
            gtilde{t + 1}(:, nz)        = statst{t}.B{nz} \ vec;
            gtilde_pre{t + 1}(:, nz)    = statst_benchmark{t}.B{nz} \ vec_pre;
        
        end
        
        gt{t + 1}           = reshape(gtilde{t + 1}, p.Nb * p.Na * p.Nz, 1);
        gt{t + 1}           = reshape(gt{t + 1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);  
        dist_vec            = reshape(gt{t+1}, p.Nb * p.Na * p.Nz, 1);
        wt{t + 1}           = p.dtildea_vec .* p.dtildeb_vec .* dist_vec;
        wt{t + 1}           = wt{t+1} / sum(wt{t+1}); 
        
        
        gt_pre{t + 1}           = reshape(gtilde_pre{t + 1}, p.Nb * p.Na * p.Nz, 1);
        gt_pre{t + 1}           = reshape(gt_pre{t + 1} ./ (p.dtildea_vec .* p.dtildeb_vec), p.Nb, p.Na, p.Nz);
        dist_pre_vec            = reshape(gt_pre{t+1}, p.Nb * p.Na * p.Nz, 1);
        wt_pre{t + 1}           = p.dtildea_vec .* p.dtildeb_vec .* dist_pre_vec;
        wt_pre{t + 1}           = wt_pre{t+1} / sum(wt_pre{t+1});  
        
        cons_vec            = reshape(statst{t+1}.cpol, p.Nb * p.Na * p.Nz, 1);
        h_vec               = reshape(statst{t+1}.hpol, p.Nb * p.Na * p.Nz, 1);
        dep_vec             = reshape(statst{t+1}.dpol, p.Nb * p.Na * p.Nz, 1);
        s_vec               = reshape(statst{t+1}.spol, p.Nb * p.Na * p.Nz, 1);
        m_vec               = reshape(statst{t+1}.mpol, p.Nb * p.Na * p.Nz, 1);
        labor_inc_vec       = reshape(statst{t+1}.labor_inc, p.Nb * p.Na * p.Nz, 1);
        liq_inc_vec         = reshape(statst{t+1}.liq_inc, p.Nb * p.Na * p.Nz, 1);
        illiq_inc_vec       = reshape(statst{t+1}.illiq_inc, p.Nb * p.Na * p.Nz, 1);
        disposable_inc_vec  = reshape(statst{t+1}.disposable_inc, p.Nb * p.Na * p.Nz, 1);

        
        cons_pre_vec            = reshape(statst_benchmark{t+1}.cpol, p.Nb * p.Na * p.Nz, 1);
        h_pre_vec               = reshape(statst_benchmark{t+1}.hpol, p.Nb * p.Na * p.Nz, 1);
        dep_pre_vec             = reshape(statst_benchmark{t+1}.dpol, p.Nb * p.Na * p.Nz, 1);
        s_pre_vec               = reshape(statst_benchmark{t+1}.spol, p.Nb * p.Na * p.Nz, 1);
        m_pre_vec               = reshape(statst_benchmark{t+1}.mpol, p.Nb * p.Na * p.Nz, 1);
        labor_inc_pre_vec       = reshape(statst_benchmark{t+1}.labor_inc, p.Nb * p.Na * p.Nz, 1);
        liq_inc_pre_vec         = reshape(statst_benchmark{t+1}.liq_inc, p.Nb * p.Na * p.Nz, 1);
        illiq_inc_pre_vec       = reshape(statst_benchmark{t+1}.illiq_inc, p.Nb * p.Na * p.Nz, 1);
        disposable_inc_pre_vec      = reshape(statst_benchmark{t+1}.disposable_inc, p.Nb * p.Na * p.Nz, 1);

        
        statst{t+1}.liq.cons_irf_percentile(k)          = (sum(wt{t + 1} .* cons_vec) - sum(wt_pre{t + 1} .* cons_pre_vec)) ./ sum(wt_pre{t + 1} .* cons_pre_vec) * 100;
        statst{t+1}.liq.h_irf_percentile(k)             = (sum(wt{t + 1} .* h_vec) - sum(wt_pre{t + 1} .* h_pre_vec)) ./ sum(wt_pre{t + 1} .* h_pre_vec) * 100;
        statst{t+1}.liq.dep_irf_percentile(k)           = (sum(wt{t + 1} .* dep_vec) - sum(wt_pre{t + 1} .* dep_pre_vec));
        statst{t+1}.liq.s_irf_percentile(k)             = (sum(wt{t + 1} .* s_vec) - sum(wt_pre{t + 1} .* s_pre_vec)) ;
        statst{t+1}.liq.m_irf_percentile(k)             = (sum(wt{t + 1} .* m_vec) - sum(wt_pre{t + 1} .* m_pre_vec)) ;
        statst{t+1}.liq.labor_inc_irf_percentile(k)     = (sum(wt{t + 1} .* labor_inc_vec) - sum(wt_pre{t + 1} .* labor_inc_pre_vec)) ./ sum(wt_pre{t + 1} .* labor_inc_pre_vec) * 100;
        statst{t+1}.liq.liq_inc_irf_percentile(k)       = (sum(wt{t + 1} .* liq_inc_vec) - sum(wt_pre{t + 1} .* liq_inc_pre_vec)) ;
        statst{t+1}.liq.illiq_inc_irf_percentile(k)     = (sum(wt{t + 1} .* illiq_inc_vec) - sum(wt_pre{t + 1} .* illiq_inc_pre_vec));
        statst{t+1}.liq.a_irf_percentile(k)             = (sum(wt{t + 1} .* a_vec) - sum(wt_pre{t + 1} .* a_vec)) ;
        statst{t+1}.liq.b_irf_percentile(k)             = (sum(wt{t + 1} .* b_vec) - sum(wt_pre{t + 1} .* b_vec));
        statst{t+1}.liq.disposable_inc_irf_percentile(k)= (sum(wt{t + 1} .* disposable_inc_vec) - sum(wt_pre{t + 1} .* disposable_inc_pre_vec)) ./sum(wt_pre{t + 1} .* disposable_inc_pre_vec) * 100;

    end
    
    end


 end