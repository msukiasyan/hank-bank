function [distr_response] = calc_distr_response(opt, glob, p,stats, statst, statst_benchmark)


    
    %% distributional IMPACT
    %% do it with respect to illiqud distribution first 
     
    dist_pre_vec    = reshape(statst_benchmark{1}.dst, p.Nb * p.Na * p.Nz, 1);   
    dist_vec        = reshape(statst{1}.dst, p.Nb * p.Na * p.Nz, 1);
    
    % percentiles at which to evaluate
    factor = statst{1}.TS /statst_benchmark{1}.TS;
    [a_percentiles_pre]  = prctilew(p.a, statst_benchmark{1}.a_marg, linspace(0.05,0.95,19));
    [b_percentiles]      = prctilew(p.b, statst_benchmark{1}.b_marg, linspace(0.05,0.95,19));
    b_percentiles        = [0,b_percentiles]
    [a_percentiles]      = factor * a_percentiles_pre;

     a_vec = p.afrombaz; 
     b_vec = p.bfrombaz; 
    % need to align the grid and the distribution closer
    
     
     cpol_pre       = griddedInterpolant({p.b,p.a,p.z}, statst_benchmark{1}.cpol);
     hpol_pre       = griddedInterpolant({p.b,p.a,p.z}, statst_benchmark{1}.hpol);
     dpol_pre       = griddedInterpolant({p.b,p.a,p.z}, statst_benchmark{1}.dpol);
     mpol_pre       = griddedInterpolant({p.b,p.a,p.z}, statst_benchmark{1}.mpol);
     spol_pre       = griddedInterpolant({p.b,p.a,p.z}, statst_benchmark{1}.spol);
     labor_inc_pre  = griddedInterpolant({p.b,p.a,p.z}, reshape(statst_benchmark{1}.labor_inc,p.Nb,p.Na,p.Nz));
     liq_inc_pre    = griddedInterpolant({p.b,p.a,p.z}, reshape(statst_benchmark{1}.liq_inc,p.Nb,p.Na,p.Nz));
     illiq_inc_pre  = griddedInterpolant({p.b,p.a,p.z}, reshape(statst_benchmark{1}.illiq_inc,p.Nb,p.Na,p.Nz));
     disposable_inc_pre  = griddedInterpolant({p.b,p.a,p.z}, reshape(statst_benchmark{1}.disposable_inc,p.Nb,p.Na,p.Nz));
     wt_pre         = griddedInterpolant({p.b,p.a,p.z}, reshape(p.dtildea_vec .* p.dtildeb_vec .* dist_pre_vec, p.Nb, p.Na, p.Nz));
     
     cpol           = griddedInterpolant({p.b,p.a,p.z}, statst{1}.cpol);
     hpol           = griddedInterpolant({p.b,p.a,p.z}, statst{1}.hpol);
     dpol           = griddedInterpolant({p.b,p.a,p.z}, statst{1}.dpol);
     mpol           = griddedInterpolant({p.b,p.a,p.z}, statst{1}.mpol);
     spol           = griddedInterpolant({p.b,p.a,p.z}, statst{1}.spol);
     labor_inc      = griddedInterpolant({p.b,p.a,p.z}, reshape(statst{1}.labor_inc,p.Nb,p.Na,p.Nz));
     liq_inc        = griddedInterpolant({p.b,p.a,p.z}, reshape(statst{1}.liq_inc,p.Nb,p.Na,p.Nz));
     illiq_inc      = griddedInterpolant({p.b,p.a,p.z}, reshape(statst{1}.illiq_inc,p.Nb,p.Na,p.Nz));
     disposable_inc  = griddedInterpolant({p.b,p.a,p.z}, reshape(statst{1}.disposable_inc,p.Nb,p.Na,p.Nz));

     wt             = griddedInterpolant({p.b,p.a,p.z}, reshape(p.dtildea_vec .* p.dtildeb_vec .* dist_vec, p.Nb, p.Na, p.Nz));
     
     cons_pre_array          = cpol_pre(({p.b,a_percentiles_pre,p.z}));
     h_pre_array             = hpol_pre(({p.b,a_percentiles_pre,p.z}));
     d_pre_array             = dpol_pre(({p.b,a_percentiles_pre,p.z}));
     m_pre_array             = mpol_pre(({p.b,a_percentiles_pre,p.z}));
     s_pre_array             = spol_pre(({p.b,a_percentiles_pre,p.z}));
     labor_inc_pre_array     = labor_inc_pre(({p.b,a_percentiles_pre,p.z}));
     liq_inc_pre_array       = liq_inc_pre(({p.b,a_percentiles_pre,p.z}));
     illiq_inc_pre_array     = illiq_inc_pre(({p.b,a_percentiles_pre,p.z}));
     disposable_inc_pre_array     = disposable_inc_pre(({p.b,a_percentiles_pre,p.z}));
     

     cons_array          = cpol(({p.b,a_percentiles,p.z}));
     h_array             = hpol(({p.b,a_percentiles,p.z}));
     d_array             = dpol(({p.b,a_percentiles,p.z}));
     m_array             = mpol(({p.b,a_percentiles,p.z}));
     s_array             = spol(({p.b,a_percentiles,p.z}));
     labor_inc_array     = labor_inc(({p.b,a_percentiles,p.z}));
     liq_inc_array       = liq_inc(({p.b,a_percentiles,p.z}));
     illiq_inc_array     = illiq_inc(({p.b,a_percentiles,p.z}));
     disposable_inc_array     = disposable_inc(({p.b,a_percentiles,p.z}));

     wt_pre_array            = wt_pre(({p.b,a_percentiles_pre,p.z})) ./ sum(sum(sum(wt_pre(({p.b,a_percentiles_pre,p.z})))));
     wt_array                = wt(({p.b,a_percentiles,p.z})) ./ sum(sum(sum(wt(({p.b,a_percentiles,p.z})))));
     
     cons_pre_vec        = sum(wt_pre_array .* cons_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);
     h_pre_vec           = sum(wt_pre_array .* h_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);
     d_pre_vec           = sum(wt_pre_array .* d_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);
     m_pre_vec           = sum(wt_pre_array .* m_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);
     s_pre_vec           = sum(wt_pre_array .* s_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);
     labor_inc_pre_vec   = sum(wt_pre_array .* labor_inc_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);
     liq_inc_pre_vec     = sum(wt_pre_array .* liq_inc_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);
     illiq_inc_pre_vec   = sum(wt_pre_array .* illiq_inc_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);
     disposable_inc_pre_vec   = sum(wt_pre_array .* disposable_inc_pre_array , [1, 3]) ./ sum(wt_pre_array, [1, 3]);

     cons_vec        = sum(wt_array .* cons_array , [1, 3]) ./ sum(wt_array, [1, 3]);
     h_vec           = sum(wt_array .* h_array , [1, 3]) ./ sum(wt_array, [1, 3]);
     d_vec           = sum(wt_array .* d_array , [1, 3]) ./ sum(wt_array, [1, 3]);
     m_vec           = sum(wt_array .* m_array , [1, 3]) ./ sum(wt_array, [1, 3]);
     s_vec           = sum(wt_array .* s_array , [1, 3]) ./ sum(wt_array, [1, 3]);
     labor_inc_vec   = sum(wt_array .* labor_inc_array , [1, 3]) ./ sum(wt_array, [1, 3]);
     liq_inc_vec     = sum(wt_array .* liq_inc_array , [1, 3]) ./ sum(wt_array, [1, 3]);
     illiq_inc_vec   = sum(wt_array .* illiq_inc_array , [1, 3]) ./ sum(wt_array, [1, 3]);
     disposable_inc_vec   = sum(wt_array .* disposable_inc_array , [1, 3]) ./ sum(wt_array, [1, 3]);

     cons_response      = (cons_vec - cons_pre_vec)./cons_pre_vec * 100;
     h_response         = (h_vec - h_pre_vec)./h_pre_vec * 100;
     d_response         = (d_vec - d_pre_vec);
     m_response         = (m_vec - m_pre_vec);
     s_response         = (s_vec - s_pre_vec);
     labor_inc_response = (labor_inc_vec - labor_inc_pre_vec)./labor_inc_pre_vec * 100;
     liq_inc_response   = (liq_inc_vec - liq_inc_pre_vec);
     illiq_inc_response = (illiq_inc_vec - illiq_inc_pre_vec);
     disposable_inc_response = (disposable_inc_vec - disposable_inc_pre_vec)./disposable_inc_pre_vec * 100;

     % now need to adjust for people with 0
     wt0                = wt_pre(({p.b,0,p.z})) ./ sum(sum(sum(wt_pre(({p.b,0,p.z})))));
     cons0_pre          = sum(wt0 .* cpol_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     h0_pre             = sum(wt0 .* hpol_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     d0_pre             = sum(wt0 .* dpol_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     m0_pre             = sum(wt0 .* mpol_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     s0_pre             = sum(wt0 .* spol_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     labor_inc0_pre     = sum(wt0 .* labor_inc_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     liq_inc0_pre       = sum(wt0 .* liq_inc_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     illiq_inc0_pre     = sum(wt0 .* illiq_inc_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     disposable_inc0_pre= sum(wt0 .* disposable_inc_pre(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);


     cons0              = sum(wt0 .* cpol(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     h0                 = sum(wt0 .* hpol(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     d0                 = sum(wt0 .* dpol(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     m0                 = sum(wt0 .* mpol(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     s0                 = sum(wt0 .* spol(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     labor_inc0         = sum(wt0 .* labor_inc(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     liq_inc0           = sum(wt0 .* liq_inc(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     illiq_inc0         = sum(wt0 .* illiq_inc(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);
     disposable_inc0    = sum(wt0 .* disposable_inc(({p.b,0,p.z})) , [1, 3]) ./ sum(wt0, [1, 3]);

     distr_response.illiq.cons          = [(cons0 - cons0_pre)./cons0_pre * 100,cons_response]; 
     distr_response.illiq.h             = [(h0 - h0_pre)./h0_pre * 100,h_response]; 
     distr_response.illiq.d             = [(d0 - d0_pre),d_response]; 
     distr_response.illiq.m             = [(m0 - m0_pre),m_response]; 
     distr_response.illiq.s             = [(s0 - s0_pre),s_response]; 
     distr_response.illiq.labor_inc      = [(labor_inc0 - labor_inc0_pre)./labor_inc0_pre * 100,labor_inc_response]; 
     distr_response.illiq.liq_inc        = [(liq_inc0 - liq_inc0_pre),liq_inc_response]; 
     distr_response.illiq.illiq_inc      = [(illiq_inc0 - illiq_inc0_pre),illiq_inc_response]; 
     distr_response.illiq.disposable_inc = [(disposable_inc0 - disposable_inc0_pre)./disposable_inc0_pre * 100,disposable_inc_response]; 

     % also: calculate distribution of liquid conditional on illiquid
     % percentiles, this will help in the interpretation
     wt_array                = wt(({p.b,[0,a_percentiles],p.z})) ./ sum(sum(sum(wt(({p.b,[0,a_percentiles],p.z})))));

          
     liquid              = griddedInterpolant({p.b,p.a,p.z}, p.bbb);
     liquid_array        = liquid(({p.b,[0,a_percentiles],p.z}));
     distr_response.illiq.liquid              = sum(wt_array .* liquid(({p.b,[0,a_percentiles],p.z})) , [1, 3]) ./ sum(wt_array, [1, 3]);
    
     %% now redo it with respect to liquid distribution
     cons_pre_array          = cpol_pre(({b_percentiles,p.a,p.z}));
     h_pre_array             = hpol_pre(({b_percentiles,p.a,p.z}));
     d_pre_array             = dpol_pre(({b_percentiles,p.a,p.z}));
     m_pre_array             = mpol_pre(({b_percentiles,p.a,p.z}));
     s_pre_array             = spol_pre(({b_percentiles,p.a,p.z}));
     labor_inc_pre_array     = labor_inc_pre(({b_percentiles,p.a,p.z}));
     liq_inc_pre_array       = liq_inc_pre(({b_percentiles,p.a,p.z}));
     illiq_inc_pre_array     = illiq_inc_pre(({b_percentiles,p.a,p.z}));
     disposable_inc_pre_array     = disposable_inc_pre(({b_percentiles,p.a,p.z}));
     cons_array          = cpol(({b_percentiles,factor*p.a,p.z}));
     h_array             = hpol(({b_percentiles,factor*p.a,p.z}));
     d_array             = dpol(({b_percentiles,factor*p.a,p.z}));
     m_array             = mpol(({b_percentiles,factor*p.a,p.z}));
     s_array             = spol(({b_percentiles,factor*p.a,p.z}));
     labor_inc_array     = labor_inc(({b_percentiles,factor*p.a,p.z}));
     liq_inc_array       = liq_inc(({b_percentiles,factor*p.a,p.z}));
     illiq_inc_array     = illiq_inc(({b_percentiles,factor*p.a,p.z}));
      disposable_inc_array     = disposable_inc(({b_percentiles,factor*p.a,p.z}));

     wt_pre_array            = wt_pre(({b_percentiles,p.a,p.z})) ./ sum(sum(sum(wt_pre(({b_percentiles,p.a,p.z})))));
     wt_array                = wt(({b_percentiles,factor*p.a,p.z})) ./ sum(sum(sum(wt(({b_percentiles,factor*p.a,p.z})))));
     
     cons_pre_vec        = sum(wt_pre_array .* cons_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);
     h_pre_vec           = sum(wt_pre_array .* h_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);
     d_pre_vec           = sum(wt_pre_array .* d_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);
     m_pre_vec           = sum(wt_pre_array .* m_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);
     s_pre_vec           = sum(wt_pre_array .* s_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);
     labor_inc_pre_vec   = sum(wt_pre_array .* labor_inc_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);
     liq_inc_pre_vec     = sum(wt_pre_array .* liq_inc_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);
     illiq_inc_pre_vec   = sum(wt_pre_array .* illiq_inc_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);
     disposable_inc_pre_vec   = sum(wt_pre_array .* disposable_inc_pre_array , [2, 3]) ./ sum(wt_pre_array, [2, 3]);

     cons_vec        = sum(wt_array .* cons_array , [2, 3]) ./ sum(wt_array, [2, 3]);
     h_vec           = sum(wt_array .* h_array , [2, 3]) ./ sum(wt_array, [2, 3]);
     d_vec           = sum(wt_array .* d_array , [2, 3]) ./ sum(wt_array, [2, 3]);
     m_vec           = sum(wt_array .* m_array , [2, 3]) ./ sum(wt_array, [2, 3]);
     s_vec           = sum(wt_array .* s_array , [2, 3]) ./ sum(wt_array, [2, 3]);
     labor_inc_vec   = sum(wt_array .* labor_inc_array , [2, 3]) ./ sum(wt_array, [2, 3]);
     liq_inc_vec     = sum(wt_array .* liq_inc_array , [2, 3]) ./ sum(wt_array, [2, 3]);
     illiq_inc_vec   = sum(wt_array .* illiq_inc_array , [2, 3]) ./ sum(wt_array, [2, 3]);
     disposable_inc_vec   = sum(wt_array .* disposable_inc_array , [2, 3]) ./ sum(wt_array, [2, 3]);

     distr_response.liq.cons      = (cons_vec - cons_pre_vec)./cons_pre_vec * 100;
     distr_response.liq.h         = (h_vec - h_pre_vec)./h_pre_vec * 100;
     distr_response.liq.d         = (d_vec - d_pre_vec);
     distr_response.liq.m         = (m_vec - m_pre_vec);
     distr_response.liq.s         = (s_vec - s_pre_vec);
     distr_response.liq.labor_inc = (labor_inc_vec - labor_inc_pre_vec)./labor_inc_pre_vec * 100;
     distr_response.liq.liq_inc   = (liq_inc_vec - liq_inc_pre_vec);
     distr_response.liq.illiq_inc = (illiq_inc_vec - illiq_inc_pre_vec);
     distr_response.liq.disposable_inc = (disposable_inc_vec - disposable_inc_pre_vec)./disposable_inc_pre_vec * 100;
 
     % also: calculate distribution of illiquid conditional on liquid
     % percentiles, this will help in the interpretation
     illiquid              = griddedInterpolant({p.b,p.a,p.z}, p.aaa);
     illiquid_array        = illiquid(({b_percentiles,p.a,p.z}));
     distr_response.liq.illiquid              = sum(wt_array .* illiquid(({b_percentiles,p.a,p.z})) , [2, 3]) ./ sum(wt_array, [2, 3]);
    
     %%
     distr.response.illiq.a_percentiles_pre = a_percentiles_pre;
     distr.response.illiq.a_percentiles = a_percentiles;
     distr.response.liq.b_percentiles = b_percentiles;
        
     %% 3-D plot
     % adjust for 0
     wt0                = wt_pre(({b_percentiles,0,p.z})); 
     cons0_pre          = sum(wt0 .* cpol_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     h0_pre             = sum(wt0 .* hpol_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     d0_pre             = sum(wt0 .* dpol_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     m0_pre             = sum(wt0 .* mpol_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     s0_pre             = sum(wt0 .* spol_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     labor_inc0_pre     = sum(wt0 .* labor_inc_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     liq_inc0_pre       = sum(wt0 .* liq_inc_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     illiq_inc0_pre     = sum(wt0 .* illiq_inc_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     disposable_inc0_pre     = sum(wt0 .* disposable_inc_pre(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);

     cons0              = sum(wt0 .* cpol(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     h0                 = sum(wt0 .* hpol(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     d0                 = sum(wt0 .* dpol(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     m0                 = sum(wt0 .* mpol(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     s0                 = sum(wt0 .* spol(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     labor_inc0         = sum(wt0 .* labor_inc(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     liq_inc0           = sum(wt0 .* liq_inc(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     illiq_inc0         = sum(wt0 .* illiq_inc(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);
     disposable_inc0         = sum(wt0 .* disposable_inc(({b_percentiles,0,p.z})) , 3) ./ sum(wt0, 3);

     cons_pre_array          = cpol_pre(({b_percentiles,a_percentiles_pre,p.z}));
     h_pre_array             = hpol_pre(({b_percentiles,a_percentiles_pre,p.z}));
     d_pre_array             = dpol_pre(({b_percentiles,a_percentiles_pre,p.z}));
     m_pre_array             = mpol_pre(({b_percentiles,a_percentiles_pre,p.z}));
     s_pre_array             = spol_pre(({b_percentiles,a_percentiles_pre,p.z}));
     labor_inc_pre_array     = labor_inc_pre(({b_percentiles,a_percentiles_pre,p.z}));
     liq_inc_pre_array       = liq_inc_pre(({b_percentiles,a_percentiles_pre,p.z}));
     illiq_inc_pre_array     = illiq_inc_pre(({b_percentiles,a_percentiles_pre,p.z}));
     disposable_inc_pre_array     = disposable_inc_pre(({b_percentiles,a_percentiles_pre,p.z}));

     cons_array          = cpol(({b_percentiles,a_percentiles,p.z}));
     h_array             = hpol(({b_percentiles,a_percentiles,p.z}));
     d_array             = dpol(({b_percentiles,a_percentiles,p.z}));
     m_array             = mpol(({b_percentiles,a_percentiles,p.z}));
     s_array             = spol(({b_percentiles,a_percentiles,p.z}));
     labor_inc_array     = labor_inc(({b_percentiles,a_percentiles,p.z}));
     liq_inc_array       = liq_inc(({b_percentiles,a_percentiles,p.z}));
     illiq_inc_array     = illiq_inc(({b_percentiles,a_percentiles,p.z}));
     disposable_inc_array     = disposable_inc(({b_percentiles,a_percentiles,p.z}));

     wt_pre_array            = wt_pre(({b_percentiles,a_percentiles_pre,p.z})) ./ sum(sum(sum(wt_pre(({b_percentiles,a_percentiles_pre,p.z})))));
     wt_array                = wt(({b_percentiles,a_percentiles,p.z})) ./ sum(sum(sum(wt(({b_percentiles,a_percentiles,p.z})))));

     cons_pre        = sum(wt_pre_array .* cons_pre_array , 3) ./ sum(wt_pre_array,  3);
     h_pre           = sum(wt_pre_array .* h_pre_array  , 3) ./ sum(wt_pre_array,  3);
     d_pre           = sum(wt_pre_array .* d_pre_array , 3) ./ sum(wt_pre_array,  3);
     m_pre           = sum(wt_pre_array .* m_pre_array , 3) ./ sum(wt_pre_array,  3);
     s_pre           = sum(wt_pre_array .* s_pre_array  , 3) ./ sum(wt_pre_array,  3);
     labor_inc_pre  = sum(wt_pre_array .* labor_inc_pre_array  , 3) ./ sum(wt_pre_array,  3);
     liq_inc_pre     = sum(wt_pre_array .* liq_inc_pre_array  , 3) ./ sum(wt_pre_array,  3);
     illiq_inc_pre   = sum(wt_pre_array .* illiq_inc_pre_array  , 3) ./ sum(wt_pre_array,  3);
     disposable_inc_pre   = sum(wt_pre_array .* disposable_inc_pre_array  , 3) ./ sum(wt_pre_array,  3);


     cons        = sum(wt_array .* cons_array , 3) ./ sum(wt_array, 3);
     h           = sum(wt_array .* h_array , 3) ./ sum(wt_array, 3);
     d           = sum(wt_array .* d_array  , 3) ./ sum(wt_array, 3);
     m           = sum(wt_array .* m_array , 3) ./ sum(wt_array, 3);
     s           = sum(wt_array .* s_array , 3) ./ sum(wt_array, 3);
     labor_inc  = sum(wt_array .* labor_inc_array , 3) ./ sum(wt_array, 3);
     liq_inc     = sum(wt_array .* liq_inc_array , 3) ./ sum(wt_array, 3);
     illiq_inc   = sum(wt_array .* illiq_inc_array , 3) ./ sum(wt_array, 3);
     disposable_inc   = sum(wt_array .* disposable_inc_array , 3) ./ sum(wt_array, 3);



     distr_response.ab.cons      = [(cons0-cons0_pre)./cons0_pre * 100,(cons - cons_pre)./cons_pre * 100];
     distr_response.ab.h         = [(h0-h0_pre)./h0_pre * 100,(h - h_pre)./h_pre * 100];
     distr_response.ab.d         = [d0-d0_pre,(d - d_pre)];
     distr_response.ab.m         = [m0-m0_pre,(m - m_pre)];
     distr_response.ab.s         = [s0-s0_pre,(s - s_pre)];
     distr_response.ab.labor_inc = [(labor_inc0 - labor_inc0_pre)./labor_inc0_pre * 100,(labor_inc - labor_inc_pre)./labor_inc_pre * 100];
     distr_response.ab.liq_inc   = [liq_inc0 - liq_inc0_pre, liq_inc - liq_inc_pre];
     distr_response.ab.illiq_inc = [illiq_inc0 - illiq_inc0_pre,(illiq_inc - illiq_inc_pre)];
     distr_response.ab.disposable_inc = [(disposable_inc0 - disposable_inc0_pre)./disposable_inc0_pre * 100,(disposable_inc - disposable_inc_pre)./disposable_inc_pre * 100];

 end