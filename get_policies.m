function sol = get_policies(opt, glob, p)
    isvalid     = true;

    %% Check
    if p.r_F - 1/p.chi1 > 0
        disp('Warning: p.r_F - 1/chi1 > 0')
    end
    
    %% Preallocation
    BU          = cell(p.Nz, 1);
    gBu         = cell(p.Nz, 1);
    gBl         = cell(p.Nz, 1);
    cFinal      = zeros(p.Nb, p.Na, p.Nz);
    dFinal      = zeros(p.Nb, p.Na, p.Nz);
    hFinal      = zeros(p.Nb, p.Na, p.Nz);
    
    %% Return at different points in state space
    Rb          = p.r_plus .* (p.bbb > 0) + p.r_minus .* (p.bbb < 0);
    Ra          = p.r_F .* ones(p.Nb, p.Na, p.Nz);
    
    %% labor supply (GHH)
    h           = ((1 - p.xi) * p.w * p.zzz / p.disutil) .^ p.varphi;
    %h           = ((1 - p.xi) * p.w / p.disutil) .^ p.varphi; % everyone provides the same h!
    lab_income  =  p.w * p.zzz .* h;
    
    %% Initial guess
    v0          = utility(((1 - p.xi) * lab_income + max((Rb) .* p.bbb, 0)),h, opt, glob, p) / p.rho; 
    v           = v0;
    V           = v;

    %% Iterate
    for n = 1:opt.maxit
        [V, dFinal, cFinal, hFinal, u, BU]  = hjb_update(opt, glob, p, v, p.Delta);
        Vchange                     = V - v;
        v                           = V;
        dist(n)                     = max(max(max(abs(Vchange))));
        if dist(n) < p.crit
            break
        end
    end
    if opt.debug_v
        if isvalid
            disp(['Value Function, Last Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
        else
            disp('Improper Transition Matrix');
        end
    end

    %% Stationary distribution
    Ai                  = cellfun(@(x) x', BU, 'UniformOutput', false);
    gmatp               = zeros(p.Nb * p.Na, p.Nz);
    gmat                = ones(p.Nb * p.Na, p.Nz);
    lmat                = eye(p.Nz, p.Nz) + p.DeltaKFE * p.la_mat_offdiag;
    
    % Set initial distribution
    ginit               = zeros(p.Nb, p.Na, p.Nz);
    for na = 1:p.Na
        ginit(10, na, 1) = p.init_dist(na);
    end
    gmat                = reshape(ginit, p.Nb * p.Na, p.Nz);
    
    for nz = 1:p.Nz
        [gBl{nz}, gBu{nz}]  = lu((1.0 - p.DeltaKFE * p.la_mat(nz, nz)) * speye(p.Nb * p.Na) - p.DeltaKFE * Ai{nz});
    end

    for itkfe = 1:opt.maxitKFE
        for nz = 1:p.Nz
            vec             = gmat * lmat(:, nz);

            gmatp(:, nz)    = gBu{nz} \ (gBl{nz} \ vec);
            
        end
        gdist           = max(max(abs(gmatp - gmat)));
        dist(itkfe)     = gdist;
        gmat            = gmatp;
        if gdist < p.critKFE
            break;
        end
    end
    if opt.debug_v
        disp(['KFE, Last Iteration ' int2str(itkfe) ', max gchange = ' num2str(gdist)]);
    end
    
    g_tilde             = reshape(gmat, p.Nb * p.Na * p.Nz, 1);
    
    g_sum               = g_tilde' * ones(p.Nb * p.Na * p.Nz, 1);
    g_tilde             = g_tilde ./ g_sum;
    
    g_stacked           = g_tilde ./ (p.dtildea_vec .* p.dtildeb_vec);
    g_stacked(g_stacked < 0)    = 0;

    g                   = reshape(g_stacked, p.Nb, p.Na, p.Nz);
    
    %% Pack
    sol.cpol            = cFinal;
    sol.hpol            = hFinal;
    sol.dpol            = dFinal;
    sol.apol            = dFinal + p.xi * lab_income + Ra .* p.aaa;
    sol.bpol            = (1 - p.xi) *  lab_income  + Rb .* p.bbb - dFinal - adjustment_cost(dFinal, p.aaa, opt, glob, p) - cFinal;
    sol.sc              = (1 - p.xi) *  lab_income  + Rb .* p.bbb - cFinal;
    sol.sd              = - dFinal - adjustment_cost(dFinal, p.aaa, opt, glob, p);
    sol.dst             = g;
    sol.gvec            = g_stacked;
    sol.isvalid         = isvalid;
    sol.u               = u;
    sol.V               = V;
    sol.BU              = BU;
    
end