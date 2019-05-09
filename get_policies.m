function sol = get_policies(opt, p)
    isvalid     = true;

    %% Solve
    rb_neg      = p.r_minus;

    if p.r_F - 1/p.chi1 > 0
        disp('Warning: p.r_F - 1/chi1 > 0')
    end

    bb          = p.b * ones(1, p.Na);
    aa          = ones(p.Nb, 1) * p.a;
    zz          = ones(p.Na, 1) * p.z;

    bbb         = zeros(p.Nb, p.Na, p.Nz); 
    aaa         = zeros(p.Nb, p.Na, p.Nz); 
    zzz         = zeros(p.Nb, p.Na, p.Nz);
    for nz = 1:p.Nz
        bbb(:, :, nz)   = bb;
        aaa(:, :, nz)   = aa;
        zzz(:, :, nz)   = p.z(nz);
    end

    Bswitch = [
        speye(p.Nb * p.Na) * p.la_mat(1, 1), speye(p.Nb * p.Na) * p.la_mat(1, 2);
        speye(p.Nb * p.Na) * p.la_mat(2, 1), speye(p.Nb * p.Na) * p.la_mat(2, 2)];

    %Preallocation
    VbF         = zeros(p.Nb, p.Na, p.Nz);
    VbB         = zeros(p.Nb, p.Na, p.Nz);
    VaF         = zeros(p.Nb, p.Na, p.Nz);
    VaB         = zeros(p.Nb, p.Na, p.Nz);
    c           = zeros(p.Nb, p.Na, p.Nz);
    updiag      = zeros(p.Nb * p.Na, p.Nz);
    lowdiag     = zeros(p.Nb * p.Na, p.Nz);
    centdiag    = zeros(p.Nb * p.Na, p.Nz);
    AAi         = cell(p.Nz, 1);
    BBi         = cell(p.Nz, 1);
    validBB     = false(p.Nb, p.Na, p.Nz);
    validFB     = false(p.Nb, p.Na, p.Nz);
    validBF     = false(p.Nb, p.Na, p.Nz);
    vcBB        = zeros(p.Nb, p.Na, p.Nz);
    vcFB        = zeros(p.Nb, p.Na, p.Nz);
    vcBF        = zeros(p.Nb, p.Na, p.Nz);
    dFinal      = zeros(p.Nb, p.Na, p.Nz);
    curvc       = zeros(p.Nb, p.Na, p.Nz);

    % return at different points in state space
    % matrices of returns
    Rb          = p.r_plus .* (bbb > 0) + p.r_minus .* (bbb < 0);
    Ra          = p.r_F .* ones(p.Nb, p.Na, p.Nz);
    
    %INITIAL GUESS
    v0          = utility(((1 - p.xi) * p.w * zzz + (Rb) .* bbb), opt, p) / p.rho;
    v           = v0;


    for n = 1:opt.maxit
        V                       = v;   
        
        % Derivatives wrt b
        %   forward difference
        VbF(1:p.Nb - 1, :, :)   = (V(2:p.Nb, :, :) - V(1:p.Nb-1, :, :)) ./ p.dbF(1:p.Nb - 1);
        %   backward difference
        VbB(2:p.Nb, :, :)       = (V(2:p.Nb, :, :) - V(1:p.Nb - 1, :, :)) ./ p.dbB(2:p.Nb);
        VbB(1, :, :)            = utility_prime((1 - p.xi) * p.w * zzz(1, :, :) + Rb(1, :, :) .* p.bmin, opt, p); %state constraint boundary condition

        % Derivatives wrt a
        %   forward difference
        VaF(:, 1:p.Na - 1, :)   = (V(:, 2:p.Na, :) - V(:, 1:p.Na - 1, :)) ./ p.daF(1:p.Na - 1);
        %   backward difference
        VaB(:, 2:p.Na, :)       = (V(:, 2:p.Na, :) - V(:, 1:p.Na - 1, :)) ./ p.daB(2:p.Na);
        
        VbF(1:p.Nb - 1, :, :)   = max(VbF(1:p.Nb - 1, :, :), 1e-8);
        VbB(2:p.Nb, :, :)       = max(VbB(2:p.Nb, :, :), 1e-8);
        VaF(:, 1:p.Na - 1, :)   = max(VaF(:, 1:p.Na - 1, :), 1e-8);
        VaB(:, 2:p.Na, :)       = max(VaB(:, 2:p.Na, :), 1e-8);

        %useful quantities
        c_B                 = VbB .^ (-1 / p.ga);
        c_F                 = VbF .^ (-1 / p.ga);
        c_0                 = (1 - p.xi) * p.w * zzz + Rb .* bbb;
        dBB                 = two_asset_kinked_FOC(VaB, VbB, aaa, opt, p);
        dFB                 = two_asset_kinked_FOC(VaB, VbF, aaa, opt, p);
        dBF                 = two_asset_kinked_FOC(VaF, VbB, aaa, opt, p);
        dFF                 = two_asset_kinked_FOC(VaF, VbF, aaa, opt, p);
        
        %% Determine the deposit components of a and b drifts
        % BB -- a backward, b backward
        validBB             = true(p.Nb, p.Na, p.Nz);
        vcBB                = VaB .* dBB - VbB .* (dBB + ...
                                two_asset_kinked_cost(dBB, aaa, opt, p));
        validBB(:, 1, :)    = 0;
        vcBB(:, 1, :)       = -1e20;
        validBB             = validBB & (dBB + two_asset_kinked_cost(dBB, aaa, opt, p) > 0) ...
                            & (dBB <= 0) & (vcBB > 0);
        % BF -- a forward, b backward
        validBF             = true(p.Nb, p.Na, p.Nz);
        vcBF                = VaF .* dBF - VbB .* (dBF + ...
                                two_asset_kinked_cost(dBF, aaa, opt, p));
        validBF(:, p.Na, :) = 0;
        validBF(1, :, :)    = 0;
        vcBF(:, p.Na, :)    = -1e20;
        vcBF(1, :, :)       = -1e20;
        validBF             = validBF & (dBF > 0) & (vcBF > 0);
        % FB -- a backward, b forward
        validFB             = true(p.Nb, p.Na, p.Nz);
        vcFB                = VaB .* dFB - VbF .* (dFB + ...
                                two_asset_kinked_cost(dFB, aaa, opt, p));
        validFB(:, 1, :)    = 0;
        vcFB(:, 1, :)       = -1e20;
        validFB(p.Nb, :, :)    = 0;
        vcFB(p.Nb, :, :)       = -1e20;
        validFB             = validFB & (dFB + two_asset_kinked_cost(dFB, aaa, opt, p) <= 0) ... 
                            & (vcFB > 0);
                        
        dFinal              = zeros(p.Nb, p.Na, p.Nz);
        curvc               = ones(p.Nb, p.Na, p.Nz) * (0);
        
        indmat              = validBB == 1;
        curvc(indmat)       = vcBB(indmat);
        dFinal(indmat)      = dBB(indmat);
        indmat              = validBF == 1 & vcBF >= curvc;
        curvc(indmat)       = vcBF(indmat);
        dFinal(indmat)      = dBF(indmat);
        indmat              = validFB == 1 & vcFB >= curvc;
        curvc(indmat)       = vcFB(indmat);
        dFinal(indmat)      = dFB(indmat);
        
        
        sdFinal             = -(dFinal + two_asset_kinked_cost(dFinal, aaa, opt, p));
        
        %% Determine consumption
        sc_B                = (1 - p.xi) * p.w * zzz + Rb .* bbb - c_B;
        validB              = true(p.Nb, p.Na, p.Nz);
        vcB                 = utility(c_B, opt, p) + VbB .* sc_B;
        validB(1, :, :)     = 0;
        vcB(1, :, :)        = -1e20;
        validB              = validB & (sc_B < 0);
        
        sc_F                = (1 - p.xi) * p.w * zzz + Rb .* bbb - c_F;
        validF              = true(p.Nb, p.Na, p.Nz);
        vcF                 = utility(c_F, opt, p) + VbF .* sc_F;
        validF(p.Nb, :, :)  = 0;
        vcF(p.Nb, :, :)     = -1e20;
        validF              = validF & (sc_F > 0);
        
        sc_0                = 0;
        vc0                 = utility(c_0, opt, p);
        
        cFinal              = zeros(p.Nb, p.Na, p.Nz);
        curvc               = ones(p.Nb, p.Na, p.Nz) * (-1e20);
        
        curvc               = vc0;
        cFinal              = c_0;
        indmat              = validB == 1 & vcB >= curvc;
        curvc(indmat)       = vcB(indmat);
        cFinal(indmat)      = c_B(indmat);
        indmat              = validF == 1 & vcF >= curvc;
        curvc(indmat)       = vcF(indmat);
        cFinal(indmat)      = c_F(indmat);

        
        scFinal             = (1 - p.xi) * p.w * zzz + Rb .* bbb - cFinal;
        
        u                   = utility(cFinal, opt, p);

        %CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
        X                   = -(min(scFinal, 0) + min(sdFinal, 0)) ./ p.dbB;
        Y                   = (min(scFinal, 0) + min(sdFinal, 0)) ./ p.dbB ...
                                - (max(scFinal, 0) + max(sdFinal, 0)) ./ p.dbF;
        Z                   = (max(scFinal, 0) + max(sdFinal, 0)) ./ p.dbF;

        if any(any(abs(X(1, :, :)) > 1e-12)) || any(any(abs(Z(p.Nb, :, :)) > 1e-12))
            disp('warning');
        end
        
        for i = 1:p.Nz
            centdiag(:, i)  = reshape(Y(:, :, i), p.Nb * p.Na, 1);
        end

        lowdiag(1:p.Nb - 1, :) = X(2:p.Nb,1,:);
        updiag(2:p.Nb, :)      = Z(1:p.Nb - 1, 1, :);
        for j = 2:p.Na
            lowdiag(1:j * p.Nb, :) = [lowdiag(1:(j - 1) * p.Nb, :); squeeze(X(2:p.Nb, j, :)); zeros(1, p.Nz)];
            updiag(1:j * p.Nb, :)  = [updiag(1:(j - 1) * p.Nb, :); zeros(1, p.Nz); squeeze(Z(1:p.Nb - 1, j, :))];
        end

        for nz = 1:p.Nz
            BBi{nz}         = spdiags(centdiag(:, nz), 0, p.Nb * p.Na, p.Nb * p.Na) ...
                + spdiags([updiag(:, nz); 0], 1, p.Nb * p.Na, p.Nb * p.Na) ...
                + spdiags([lowdiag(:, nz); 0], -1, p.Nb * p.Na, p.Nb * p.Na);
        end

        % CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
        MB                  = min(dFinal, 0);
        MF                  = max(dFinal, 0) + p.xi * p.w * zzz + Ra .* aaa;
        MF(:, p.Na, :)      = 0.0;
        chi                 = -MB ./ p.daB;
        yy                  = MB ./ p.daB - MF ./ p.daF;
        zeta                = MF ./ p.daF;

        if any(any(abs(MB(:, 1, :)) > 1e-12)) || any(any(abs(MF(:, p.Na, :)) > 1e-12))
            disp('warning');
        end

        % Matrix AAi
        for nz = 1:p.Nz
            AAupdiag        = [zeros(p.Nb, 1); reshape(zeta(:, 1:p.Na - 1, nz), (p.Na - 1) * p.Nb, 1)];
            AAcentdiag      = reshape(yy(:, :, nz), p.Na * p.Nb, 1);
            AAlowdiag       = reshape(chi(:, 2:p.Na, nz), (p.Na - 1) * p.Nb, 1);

            % Add up the upper, center, and lower diagonal into a sparse matrix
            AAi{nz}         = spdiags(AAcentdiag, 0, p.Nb * p.Na, p.Nb * p.Na) ...
                + spdiags(AAlowdiag, -p.Nb, p.Nb * p.Na, p.Nb * p.Na) ...
                + spdiags(AAupdiag, p.Nb, p.Nb * p.Na, p.Nb * p.Na);

        end
        
        Ai                  = cellfun(@plus, AAi, BBi, 'UniformOutput', false);
%         if any(cellfun(@(x) full(max(abs(sum(x, 2)))), Ai) > 1e-6)
%            isvalid  = false;
%            break
%         end
        
        u_ba_z              = reshape(u, p.Nb * p.Na, p.Nz);
        V_ba_z              = reshape(V, p.Nb * p.Na, p.Nz);
        Vp_ba_z             = zeros(p.Nb * p.Na, p.Nz);
        
        for nz = 1:p.Nz
            B               = (1 + p.Delta * p.rho - p.Delta * p.la_mat(nz, nz)) * speye(p.Nb * p.Na) - p.Delta * Ai{nz};
            vec             = p.Delta * u_ba_z(:, nz) + V_ba_z(:, nz) + p.Delta * V_ba_z * p.la_mat_offdiag(nz, :)';

            Vp_ba_z(:, nz)  = B \ vec;
        end
        V                   = reshape(Vp_ba_z, p.Nb, p.Na, p.Nz);


        Vchange             = V - v;
        v                   = V;

        dist(n)             = max(max(max(abs(Vchange))));
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STATIONARY DISTRIBUTION %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    d               = dFinal;
    d(p.Nb, :, :)   = d(p.Nb - 1, :, :);
    m               = d + p.xi * p.w * zzz + Ra .* aaa;
    s               = (1 - p.xi) * p.w * zzz + Rb .* bbb - d - two_asset_kinked_cost(d, aaa, opt, p) - cFinal;
    
    m(:, 1, :)      = max(m(:, 1, :), 0.0);
    m(:, p.Na, :)   = min(m(:, p.Na, :), 0.0);
    s(1, :, :)      = max(s(1, :, :), 0.0);
    s(p.Nb, :, :)   = min(s(p.Nb, :, :), 0.0);
    
    % RECONSTRUCT TRANSITION MATRIX WITH SIMPLER UPWIND SCHEME
    X                   = -min(s, 0) ./ p.dbB;
    Y                   = min(s, 0) ./ p.dbB - max(s, 0) ./ p.dbF;
    Z                   = max(s, 0) ./ p.dbF;

    for i = 1:p.Nz
        centdiag(:, i)  = reshape(Y(:, :, i), p.Nb * p.Na, 1);
    end

    lowdiag(1:p.Nb - 1, :) = X(2:p.Nb, 1, :);
    updiag(2:p.Nb, :)      = Z(1:p.Nb - 1, 1, :);
    for j = 2:p.Na
        lowdiag(1:j * p.Nb, :) = [lowdiag(1:(j - 1) * p.Nb, :); squeeze(X(2:p.Nb, j, :)); zeros(1, p.Nz)];
        updiag(1:j * p.Nb, :)  = [updiag(1:(j - 1) * p.Nb,:); zeros(1, p.Nz); squeeze(Z(1:p.Nb - 1, j, :))];
    end

    for nz = 1:p.Nz
        BBi{nz} = spdiags(centdiag(:, nz), 0, p.Nb * p.Na, p.Nb * p.Na) ...
            + spdiags([updiag(:, nz); 0], 1, p.Nb * p.Na, p.Nb * p.Na) ...
            + spdiags([lowdiag(:, nz); 0], -1, p.Nb * p.Na , p.Nb * p.Na);
    end

    BB                  = [BBi{1}, sparse(p.Nb * p.Na, p.Nb * p.Na); sparse(p.Nb * p.Na, p.Nb * p.Na), BBi{2}];

    %CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
    chi                 = -min(m,0) ./ p.daB;
    yy                  =  min(m,0) ./ p.daB - max(m,0) ./ p.daF;
    zeta                = max(m,0) ./ p.daF;


    %MATRIX AAi
    for nz = 1:p.Nz

        AAupdiag        = [zeros(p.Nb, 1); reshape(zeta(:, 1:p.Na - 1, nz), (p.Na - 1) * p.Nb, 1)];
        AAcentdiag      = reshape(yy(:, :, nz), p.Na * p.Nb, 1);
        AAlowdiag       = reshape(chi(:, 2:p.Na, nz), (p.Na - 1) * p.Nb, 1);

        %Add up the upper, center, and lower diagonal into a sparse matrix
        AAi{nz}         = spdiags(AAcentdiag, 0, p.Nb * p.Na, p.Nb * p.Na) ...
            + spdiags(AAlowdiag, -p.Nb, p.Nb * p.Na, p.Nb * p.Na) ...
            + spdiags(AAupdiag, p.Nb, p.Nb * p.Na, p.Nb * p.Na);

    end
    
    Ai                  = cellfun(@(x, y) (x + y)', AAi, BBi, 'UniformOutput', false);
    gmatp               = zeros(p.Nb * p.Na, p.Nz);
    gmat                = ones(p.Nb * p.Na, p.Nz);
    lmat                = eye(nz, nz) + p.DeltaKFE * p.la_mat_offdiag;
    
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
    sol.dpol            = dFinal;
    sol.apol            = dFinal + p.xi * p.w * zzz + Ra .* aaa;
    sol.bpol            = (1 - p.xi) * p.w * zzz + Rb .* bbb - dFinal - two_asset_kinked_cost(dFinal, aaa, opt, p) - cFinal;
    sol.sc              = (1 - p.xi) * p.w * zzz + Rb .* bbb - cFinal;
    sol.sd              = - dFinal - two_asset_kinked_cost(dFinal, aaa, opt, p);
    sol.dst             = g;
    sol.isvalid         = isvalid;
    sol.V               = V;
    
end