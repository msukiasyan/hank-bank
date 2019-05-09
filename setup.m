function pars   = setup(opt, p)

    %% Income transition matrix
    p.la_mat_diag       = spdiags(spdiags(p.la_mat, 0), 0, p.Nz, p.Nz);
    p.la_mat_offdiag    = p.la_mat - p.la_mat_diag;

    %% Create uneven grids
    p.b                 = p.bmin + linspace(0, 1, p.Nb)' .^ p.bcurve * (p.bmax - p.bmin);
    p.a                 = p.amin + linspace(0, 1, p.Na) .^ p.acurve * (p.amax - p.amin);
    
    % Correction for a grid
    p.a(1:10)           = linspace(0, p.a(10), 10);
    
    %% Create forward and backward differences
    p.dbF               = zeros(p.Nb, 1);
    p.dbB               = zeros(p.Nb, 1);
    p.daF               = zeros(1, p.Na);
    p.daB               = zeros(1, p.Na);
    
    p.db                = (p.bmax - p.bmin) / (p.Nb - 1);
    p.dbF(1:p.Nb - 1)   = (p.b(2:p.Nb) - p.b(1:p.Nb - 1));
    p.dbF(p.Nb)         = p.dbF(p.Nb - 1);
    p.dbB(2:p.Nb)       = p.dbF(1:p.Nb - 1);
    p.dbB(1)            = p.dbB(2);

    p.da                = (p.amax-p.amin) / (p.Na - 1);
    p.daF(1:p.Na - 1)   = (p.a(2:p.Na) - p.a(1:p.Na - 1));
    p.daF(p.Na)         = p.daF(p.Na - 1);
    p.daB(2:p.Na)       = p.daF(1:p.Na - 1);
    p.daB(1)            = p.daB(2);
    
    p.dbF_baz           = zeros(p.Nb, p.Na, p.Nz);
    p.dbB_baz           = zeros(p.Nb, p.Na, p.Nz);
    p.daF_baz           = zeros(p.Nb, p.Na, p.Nz);
    p.daB_baz           = zeros(p.Nb, p.Na, p.Nz);
    
    p.dbF_ba            = p.dbF * ones(1, p.Na);
    p.dbB_ba            = p.dbB * ones(1, p.Na);
    p.daF_ba            = ones(p.Nb, 1) * p.daF;
    p.daB_ba            = ones(p.Nb, 1) * p.daB;
    
    for nz = 1:p.Nz
        p.dbF_baz(:, :, nz)     = p.dbF_ba;
        p.dbB_baz(:, :, nz)     = p.dbB_ba;
        p.daF_baz(:, :, nz)     = p.daF_ba;
        p.daB_baz(:, :, nz)     = p.daB_ba;
    end
    
    p.bfrombaz          = p.b(ind2subind([p.Nb, p.Na, p.Nz], [1:p.Nb * p.Na * p.Nz], 1));
    p.afrombaz          = p.a(ind2subind([p.Nb, p.Na, p.Nz], [1:p.Nb * p.Na * p.Nz], 2))';
    p.zfrombaz          = p.z(ind2subind([p.Nb, p.Na, p.Nz], [1:p.Nb * p.Na * p.Nz], 3))';
    
    p.bindfrombaz       = ind2subind([p.Nb, p.Na, p.Nz], [1:p.Nb * p.Na * p.Nz], 1);
    p.aindfrombaz       = ind2subind([p.Nb, p.Na, p.Nz], [1:p.Nb * p.Na * p.Nz], 2)';
    p.zindfrombaz       = ind2subind([p.Nb, p.Na, p.Nz], [1:p.Nb * p.Na * p.Nz], 3)';
    
    % Construct adjustment vectors used for non-uniform grid KFE
    p.dtildea           = zeros(1, p.Na);
    p.dtildeb           = zeros(p.Nb, 1);
    
    p.dtildea(2:p.Na)   = 0.5 * (p.daF(2:p.Na) + p.daB(2:p.Na));
    p.dtildea(1)        = 0.5 * p.daF(1);
    p.dtildea(p.Na)     = 0.5 * p.daB(p.Na);
    
    p.dtildeb(2:p.Nb)   = 0.5 * (p.dbF(2:p.Nb) + p.dbB(2:p.Nb));
    p.dtildeb(1)        = 0.5 * p.dbF(1);
    p.dtildeb(p.Nb)     = 0.5 * p.dbB(p.Nb);
    
    p.dtildea_baz       = zeros(p.Nb, p.Na, p.Nz);
    p.dtildeb_baz       = zeros(p.Nb, p.Na, p.Nz);
    
    p.dtildea_ba        = ones(p.Nb, 1) * p.dtildea;
    p.dtildeb_ba        = p.dtildeb * ones(1, p.Na);
    
    for nz = 1:p.Nz
        p.dtildea_baz(:, :, nz)     = p.dtildea_ba;
        p.dtildeb_baz(:, :, nz)     = p.dtildeb_ba;
    end
    
    p.dtildea_vec       = reshape(p.dtildea_baz, p.Nb * p.Na * p.Nz, 1);
    p.dtildeb_vec       = reshape(p.dtildeb_baz, p.Nb * p.Na * p.Nz, 1);
    
    pars                = p;
end