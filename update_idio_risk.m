function p  = update_idio_risk(opt, glob, p)
    %% Finalize shock grid
    p.z             = p.zfactor * (p.zbase - mean(p.zbase)) + mean(p.zbase);
    
    %% Find the stationary distribution for z
    [E, V]          = eig(p.la_mat);
    p.z_dist        = E(:, abs(diag(V)) < 1e-12);
    p.z_dist        = p.z_dist' / sum(p.z_dist);
    
    %% Expand the grids
    bb          = p.b * ones(1, p.Na);
    aa          = ones(p.Nb, 1) * p.a;
    zz          = ones(p.Na, 1) * p.z;

    p.bbb       = zeros(p.Nb, p.Na, p.Nz); 
    p.aaa       = zeros(p.Nb, p.Na, p.Nz); 
    p.zzz       = zeros(p.Nb, p.Na, p.Nz);
    for nz = 1:p.Nz
        p.bbb(:, :, nz)     = bb;
        p.aaa(:, :, nz)     = aa;
        p.zzz(:, :, nz)     = p.z(nz);
    end

end