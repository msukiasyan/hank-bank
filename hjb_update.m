function [V1, dFinal, cFinal, h, u, BU, B] = hjb_update(opt, glob, p, V, Delta)
    isvalid     = true;

    %% Prepare objects
    VbF         = zeros(p.Nb, p.Na, p.Nz);
    VbB         = zeros(p.Nb, p.Na, p.Nz);
    VaF         = zeros(p.Nb, p.Na, p.Nz);
    VaB         = zeros(p.Nb, p.Na, p.Nz);
    c           = zeros(p.Nb, p.Na, p.Nz);
    h           = zeros(p.Nb, p.Na, p.Nz);
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
    curvc       = zeros(p.Nb, p.Na, p.Nz);
    B           = cell(p.Nz, 1);
    BU          = cell(p.Nz, 1);

    % return at different points in state space
    % matrices of returns
    Rb          = p.r_plus .* (p.bbb > 0) + p.r_minus .* (p.bbb < 0);
    Ra          = p.r_F .* ones(p.Nb, p.Na, p.Nz);
    
    adriftb     = zeros(p.Nb, p.Na, p.Nz);
    if opt.divtoliq
        adriftb = Ra .* p.aaa;
    end
    
    
    % labor supply - easy with GHH
    
    h                   = ((1 - p.xi) * p.w * p.zzz / p.disutil) .^ p.frisch;
    lab_income          = p.w * p.zzz .* h;
    
    %% Solve
    % Derivatives wrt b
    %   forward difference
    VbF(1:p.Nb - 1, :, :)   = (V(2:p.Nb, :, :) - V(1:p.Nb-1, :, :)) ./ p.dbF(1:p.Nb - 1);
    %   backward difference
    VbB(2:p.Nb, :, :)       = (V(2:p.Nb, :, :) - V(1:p.Nb - 1, :, :)) ./ p.dbB(2:p.Nb);
    VbB(1, :, :)            = utility_prime((1 - p.xi) * lab_income(1, :, :) + Rb(1, :, :) .* p.bmin + adriftb(1, :, :), h(1, :, :), opt, glob, p); %state constraint boundary condition

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
    
    
    % first need to get labor supply
    % here it is the same for all
    % consumption levels!
    
    
    % i am a bit confused about p.xi - this is something to figure out
    % later
    
   
    
    c_B                 = VbB .^ (-1 / p.ga) + p.disutil * 1 / (1 + 1 / p.frisch) * h .^ (1 + 1 / p.frisch);
    c_F                 = VbF .^ (-1 / p.ga) + p.disutil * 1 / (1 + 1 / p.frisch) * h .^ (1 + 1 / p.frisch);
    

    

    c_0                 = (1 - p.xi) * lab_income + Rb .* p.bbb + adriftb;
    


    
    dBB                 = optimal_deposit(VaB, VbB, p.aaa, opt, glob, p);
    dFB                 = optimal_deposit(VaB, VbF, p.aaa, opt, glob, p);
    dBF                 = optimal_deposit(VaF, VbB, p.aaa, opt, glob, p);

    %% Determine the deposit components of a and b drifts
    % BB -- a backward, b backward
    validBB             = true(p.Nb, p.Na, p.Nz);
    vcBB                = VaB .* dBB - VbB .* (dBB + ...
                            adjustment_cost(dBB, p.aaa, opt, glob, p));
    validBB(:, 1, :)    = 0;
    vcBB(:, 1, :)       = -1e20;
    validBB             = validBB & (dBB + adjustment_cost(dBB, p.aaa, opt, glob, p) > 0) ...
                        & (dBB <= 0) & (vcBB > 0);
    % BF -- a forward, b backward
    validBF             = true(p.Nb, p.Na, p.Nz);
    vcBF                = VaF .* dBF - VbB .* (dBF + ...
                            adjustment_cost(dBF, p.aaa, opt, glob, p));
    validBF(:, p.Na, :) = 0;
    validBF(1, :, :)    = 0;
    vcBF(:, p.Na, :)    = -1e20;
    vcBF(1, :, :)       = -1e20;
    validBF             = validBF & (dBF > 0) & (vcBF > 0);
    % FB -- a backward, b forward
    validFB             = true(p.Nb, p.Na, p.Nz);
    vcFB                = VaB .* dFB - VbF .* (dFB + ...
                            adjustment_cost(dFB, p.aaa, opt, glob, p));
    validFB(:, 1, :)    = 0;
    vcFB(:, 1, :)       = -1e20;
    validFB(p.Nb, :, :)    = 0;
    vcFB(p.Nb, :, :)       = -1e20;
    validFB             = validFB & (dFB + adjustment_cost(dFB, p.aaa, opt, glob, p) <= 0) ... 
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


    sdFinal             = -(dFinal + adjustment_cost(dFinal, p.aaa, opt, glob, p));

    %% Determine consumption
    sc_B                = (1 - p.xi) * lab_income + Rb .* p.bbb - c_B + adriftb;
    validB              = true(p.Nb, p.Na, p.Nz);
    vcB                 = utility(c_B, h, opt, glob, p) + VbB .* sc_B;
    validB(1, :, :)     = 0;
    vcB(1, :, :)        = -1e20;
    validB              = validB & (sc_B < 0);

    sc_F                = (1 - p.xi) * lab_income + Rb .* p.bbb - c_F + adriftb;
    validF              = true(p.Nb, p.Na, p.Nz);
    vcF                 = utility(c_F, h, opt, glob, p) + VbF .* sc_F;
    validF(p.Nb, :, :)  = 0;
    vcF(p.Nb, :, :)     = -1e20;
    validF              = validF & (sc_F > 0);

    sc_0                = 0;
    vc0                 = utility(c_0, h, opt, glob, p);

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


    scFinal             = (1 - p.xi) * lab_income + Rb .* p.bbb - cFinal + adriftb;

    u                   = utility(cFinal, h, opt, glob, p);

    %CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
    X                   = -(min(scFinal, 0) + min(sdFinal, 0)) ./ p.dbB;
    Y                   = (min(scFinal, 0) + min(sdFinal, 0)) ./ p.dbB ...
                            - (max(scFinal, 0) + max(sdFinal, 0)) ./ p.dbF;
    Z                   = (max(scFinal, 0) + max(sdFinal, 0)) ./ p.dbF;

    if any(any(abs(X(1, :, :)) > 1e-12)) || any(any(abs(Z(p.Nb, :, :)) > 1e-12))
        % disp('warning');
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
    MF                  = max(dFinal, 0) + p.xi * lab_income;
    if ~opt.divtoliq
        MF              = MF + Ra .* p.aaa;
    end
    MF(:, p.Na, :)      = 0.0;
    
    chi                 = -MB ./ p.daB;
    yy                  = MB ./ p.daB - MF ./ p.daF;
    zeta                = MF ./ p.daF;

    if any(any(abs(MB(:, 1, :)) > 1e-12)) || any(any(abs(MF(:, p.Na, :)) > 1e-12))
        % disp('warning');
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
    u_ba_z              = reshape(u, p.Nb * p.Na, p.Nz);
    V_ba_z              = reshape(V, p.Nb * p.Na, p.Nz);
    Vp_ba_z             = zeros(p.Nb * p.Na, p.Nz);

    for nz = 1:p.Nz
        B{nz}           = (1 + Delta * p.rho - Delta * p.la_mat(nz, nz)) * speye(p.Nb * p.Na) - Delta * Ai{nz};
        vec             = Delta * u_ba_z(:, nz) + V_ba_z(:, nz) + Delta * V_ba_z * p.la_mat_offdiag(nz, :)';

        Vp_ba_z(:, nz)  = B{nz} \ vec;
    end
    V1                  = reshape(Vp_ba_z, p.Nb, p.Na, p.Nz);
    
    %% If the transition matrix is requested
    if nargout > 4                                                          
        d               = dFinal;
        d(p.Nb, :, :)   = d(p.Nb - 1, :, :);
        m               = d + p.xi * lab_income;
        s               = (1 - p.xi) * lab_income + Rb .* p.bbb - d - adjustment_cost(d, p.aaa, opt, glob, p) - cFinal;

        if ~opt.divtoliq
            m           = m + Ra .* p.aaa;
        else
            s           = s + Ra .* p.aaa;
        end
        
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
        BU                  = cellfun(@(x, y) (x + y), AAi, BBi, 'UniformOutput', false);
    end
end