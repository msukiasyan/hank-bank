function sol = get_policies(opt, p)
    isvalid     = true;
    Delta       = p.Delta;
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

    %return at different points in state space
    %matrix of liquid returns
    Rb          = p.r_plus .* (bbb > 0) + p.r_minus .* (bbb < 0);
    raa         = p.r_F .* ones(1, p.Na);
    %if p.r_F>>rb, impose tax on p.r_F*a at high a, otherwise some households
    %accumulate infinite illiquid wealth (not needed if p.r_F is close to or less than rb)
    tau         = 10; 
    raa         = p.r_F .* (1 - (1.33 .* p.amax ./ p.a) .^ (1 - tau));
    %matrix of illiquid returns
    Ra(:, :, 1) = ones(p.Nb, 1) * raa;
    Ra(:, :, 2) = ones(p.Nb, 1) * raa;
    
    %INITIAL GUESS
    %v0 = (((1-p.xi)*p.w*zzz + p.r_F.*aaa + (Rb) .*bbb).^(1-p.ga))/(1-p.ga)/p.rho;
    v0          = (((1 - p.xi) * p.w * zzz + (Rb) .* bbb) .^ (1 - p.ga)) / (1 - p.ga) / p.rho;
    v           = v0;


    for n = 1:opt.maxit
        Delta                   = p.Delta;
        
        V                       = v;   
        
        %DERIVATIVES W.R.T. b
        % forward difference
        VbF(1:p.Nb - 1, :, :)   = (V(2:p.Nb, :, :) - V(1:p.Nb-1, :, :)) ./ p.dbF(1:p.Nb - 1);
        VbF(p.Nb, :, :)         = ((1 - p.xi) * p.w * zzz(p.Nb, :, :) + Rb(p.Nb, :, :) .* p.bmax) .^ (-p.ga); %state constraint boundary condition
        % backward difference
        VbB(2:p.Nb, :, :)       = (V(2:p.Nb, :, :) - V(1:p.Nb - 1, :, :)) ./ p.dbB(2:p.Nb);
        VbB(1, :, :)            = ((1 - p.xi) * p.w * zzz(1, :, :) + Rb(1, :, :) .* p.bmin) .^ (-p.ga); 
                                    %state constraint boundary condition

        %DERIVATIVES W.R.T. a
        % forward difference
        VaF(:, 1:p.Na - 1, :)  = (V(:, 2:p.Na, :) - V(:, 1:p.Na - 1, :)) ./ p.daF(1:p.Na - 1);
        % backward difference
        VaB(:, 2:p.Na, :)      = (V(:, 2:p.Na, :) - V(:, 1:p.Na - 1, :)) ./ p.daB(1:p.Na - 1);

        %useful quantities
        c_B                 = max(VbB, 1e-7) .^ (-1 / p.ga);
        c_F                 = max(VbF, 1e-7) .^ (-1 / p.ga);
        dBB                 = two_asset_kinked_FOC(VaB, VbB, aaa, opt, p);
        dFB                 = two_asset_kinked_FOC(VaB, VbF, aaa, opt, p);
        %VaF(:,p.Na,:) = VbB(:,p.Na,:).*(1-p.r_F.*chi1 - chi1*p.w*zzz(:,p.Na,:)./a(:,p.Na,:));
        dBF                 = two_asset_kinked_FOC(VaF, VbB, aaa, opt, p);
        %VaF(:,p.Na,:) = VbF(:,p.Na,:).*(1-p.r_F.*chi1 - chi1*p.w*zzz(:,p.Na,:)./a(:,p.Na,:));
        dFF                 = two_asset_kinked_FOC(VaF, VbF, aaa, opt, p);
        
        %% Determine the deposit components of a and b drifts
        % BB -- a backward, b backward
        validBB             = true(p.Nb, p.Na, p.Nz);
        vcBB                = VaB .* dBB - VbB .* (dBB + ...
                                two_asset_kinked_cost(dBB, aaa, opt, p));
        validBB(:, 1, :)    = 0;
        vcBB(:, 1, :)       = -1e12;
        validBB             = validBB & (dBB + two_asset_kinked_cost(dBB, aaa, opt, p) > 0) ...
                            & (dBB <= 0) & (vcBB > 0);
        % BF -- a forward, b backward
        validBF             = true(p.Nb, p.Na, p.Nz);
        vcBF                = VaF .* dBF - VbB .* (dBF + ...
                                two_asset_kinked_cost(dBF, aaa, opt, p));
        validBF(:, p.Na, :) = 0;
        validBF(1, :, :)    = 0;
        vcBF(:, p.Na, :)    = -1e12;
        vcBF(1, :, :)       = -1e12;
        validBF             = validBF & (dBF > 0) & (vcBF > 0);
        % FB -- a backward, b forward
        validFB             = true(p.Nb, p.Na, p.Nz);
        vcFB                = VaB .* dFB - VbF .* (dFB + ...
                                two_asset_kinked_cost(dFB, aaa, opt, p));
        validFB(:, 1, :)    = 0;
        vcFB(:, 1, :)       = -1e12;
        validFB(p.Nb, :, :)    = 0;
        vcFB(p.Nb, :, :)       = -1e12;
        validFB             = validFB & (dFB + two_asset_kinked_cost(dFB, aaa, opt, p) <= 0) ... 
                            & (vcFB > 0);
                        
        dFinal              = zeros(p.Nb, p.Na, p.Nz);
        curvc               = zeros(p.Nb, p.Na, p.Nz);
        
        indmat              = validBB == 1;
        curvc(indmat)       = vcBB(indmat);
        dFinal(indmat)      = dBB(indmat);
        indmat              = validBF == 1 & vcBF > curvc;
        curvc(indmat)       = vcBF(indmat);
        dFinal(indmat)      = dBF(indmat);
        indmat              = validFB == 1 & vcFB > curvc;
        curvc(indmat)       = vcFB(indmat);
        dFinal(indmat)      = dFB(indmat);
        
        sdFinal             = -(dFinal + two_asset_kinked_cost(dFinal, aaa, opt, p));
        
        %% 
        %UPWIND SCHEME
        d_B                 = (dBF > 0) .* dBF + (dBB < 0) .* dBB;
        %state constraints at amin and p.amax
        d_B(:, 1, :)        = (dBF(:, 1, :) > 1e-12) .* dBF(:, 1, :); %make sure d>=0 at p.amax, don't use VaB(:,1,:)
        d_B(:, p.Na, :)     = (dBB(:, p.Na, :) < -1e-12) .* dBB(:, p.Na, :); %make sure d<=0 at p.amax, don't use VaF(:,p.Na,:)
        d_B(1, 1, :)        = max(d_B(1, 1, :), 0);
        %split drift of b and upwind separately
        sc_B                = (1 - p.xi) * p.w * zzz + Rb .* bbb - c_B;
        sd_B                = (-d_B - two_asset_kinked_cost(d_B, aaa, opt, p));

        d_F                 = (dFF > 0) .* dFF + (dFB < 0) .* dFB;
        %state constraints at amin and p.amax
        d_F(:, 1, :)        = (dFF(:, 1, :) > 1e-12) .* dFF(:, 1, :); %make sure d>=0 at amin, don't use VaB(:,1,:)
        d_F(:, p.Na, :)        = (dFB(:, p.Na, :) < -1e-12) .* dFB(:, p.Na, :); %make sure d<=0 at p.amax, don't use VaF(:,p.Na,:)

        %split drift of b and upwind separately
        sc_F                = (1 - p.xi) * p.w * zzz + Rb .* bbb - c_F;
        sd_F                = (-d_F - two_asset_kinked_cost(d_F, aaa, opt, p));
        sd_F(p.Nb, :, :)       = min(sd_F(p.Nb, :, :), 0);

        Ic_B                = (sc_B < -1e-12);
        Ic_F                = (sc_F > 1e-12) .* (1 - Ic_B);
        Ic_0                = 1 - Ic_F - Ic_B;

        Id_F                = (sd_F > 1e-12);
        Id_B                = (sd_B < -1e-12) .* (1 - Id_F);
        Id_B(1, :, :)       = 0;
        Id_F(p.Nb, :, :)       = 0; 
        Id_B(p.Nb, :, :)       = 1; %don't use VbF at p.bmax so as not to pick up articial state constraint
        Id_0                = 1 - Id_F - Id_B;

        c_0                 = (1 - p.xi) * p.w * zzz + Rb .* bbb;

        c                   = c_F .* Ic_F + c_B .* Ic_B + c_0 .* Ic_0;
        u                   = c .^ (1 - p.ga) / (1 - p.ga);

        %CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
        X                   = -Ic_B .* sc_B ./ p.dbB - min(sdFinal, 0) ./ p.dbB;
        Y                   = (Ic_B .* sc_B + min(sdFinal, 0)) ./ p.dbB ...
                                - (Ic_F .* sc_F + max(sdFinal, 0)) ./ p.dbF;
        Z                   = Ic_F .* sc_F ./ p.dbF + max(sdFinal, 0) ./ p.dbF;

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

        BB                  = [BBi{1}, sparse(p.Nb * p.Na, p.Nb * p.Na); sparse(p.Nb * p.Na, p.Nb * p.Na), BBi{2}];

        %CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
        % dB                  = Id_B .* dBB + Id_F .* dFB;
        % dF                  = Id_B .* dBF + Id_F .* dFF;
        % d_final             = Id_B .* d_B + Id_F .* d_F;            %!!
        MB                  = min(dFinal, 0);   %!!
        MF                  = max(dFinal, 0) + p.xi * p.w * zzz + Ra .* aaa;   %!!
        MB(:, p.Na, :)      = p.xi * p.w * zzz(:, p.Na, :) + dFinal(:, p.Na, :) + Ra(:, p.Na, :) .* p.amax; %this is hopefully negative
        MF(:, p.Na, :)      = 0;
        chi                 = -MB ./ p.daB;
        yy                  = MB ./ p.daB - MF ./ p.daF;
        zeta                = MF ./ p.daF;

        %MATRIX AAi
        for nz = 1:p.Nz
            %This will be the upperdiagonal of the matrix AAi
            AAupdiag        = [zeros(p.Nb, 1); reshape(zeta(:, 1:p.Na - 1, nz), (p.Na - 1) * p.Nb, 1)];
            % AAupdiag        = zeros(p.Nb, 1); %This is necessary because of the peculiar way spdiags is defined.
            % for j=1:p.Na - 1 % MISTAKE CORRECTED: p.Na -> p.Na-1
            %     AAupdiag    = [AAupdiag; zeta(:, j, nz)];
            % end

            %This will be the center diagonal of the matrix AAi
            AAcentdiag      = reshape(yy(:, :, nz), p.Na * p.Nb, 1);
            % AAcentdiag      = yy(:, 1, nz);
            % for j = 2:p.Na - 1
            %     AAcentdiag  = [AAcentdiag; yy(:, j, nz)];
            % end
            % AAcentdiag      = [AAcentdiag; yy(:, p.Na, nz)];

            %This will be the lower diagonal of the matrix AAi
            AAlowdiag       = reshape(chi(:, 2:p.Na, nz), (p.Na - 1) * p.Nb, 1);
            % AAlowdiag       = chi(:, 2, nz);
            % for j = 3:p.Na
            %     AAlowdiag   = [AAlowdiag; chi(:, j, nz)];
            % end

            %Add up the upper, center, and lower diagonal into a sparse matrix
            AAi{nz}         = spdiags(AAcentdiag, 0, p.Nb * p.Na, p.Nb * p.Na) ...
                + spdiags(AAlowdiag, -p.Nb, p.Nb * p.Na, p.Nb * p.Na) ...
                + spdiags(AAupdiag, p.Nb, p.Nb * p.Na, p.Nb * p.Na);

        end

        AA                  = [AAi{1}, sparse(p.Nb * p.Na, p.Nb * p.Na); sparse(p.Nb * p.Na, p.Nb * p.Na), AAi{2}];

        A                   = AA + BB + Bswitch;

        if max(abs(sum(A, 2))) > 1e-7
           isvalid  = false;
           break
        end

        B                   = (1 / Delta + p.rho) * speye(p.Nb * p.Na * p.Nz) - A;

        u_stacked           = reshape(u, p.Nb * p.Na * p.Nz,1);
        V_stacked           = reshape(V, p.Nb * p.Na * p.Nz,1);

        vec                 = u_stacked + V_stacked / Delta;

        V_stacked           = B \ vec; %SOLVE SYSTEM OF EQUATIONS
        
        V                   = reshape(V_stacked, p.Nb, p.Na, p.Nz);   


        Vchange             = V - v;
        v                   = V;

        b_dist(:, n)        = max(abs(Vchange(:, :, 2)), [], 2);
        a_dist(:, n)        = max(max(abs(Vchange), [], 3), [], 1);
        ab_dist(:, :, n)    = max(abs(Vchange), [], 3);

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

    d               = dFinal;
    m               = d + p.xi * p.w * zzz + Ra .* aaa;
    s               = (1 - p.xi) * p.w * zzz + Rb .* bbb - d - two_asset_kinked_cost(d, aaa, opt, p) - c;

    sc              = (1 - p.xi) * p.w * zzz + Rb .* bbb - c;
    sd              = - d - two_asset_kinked_cost(d, aaa, opt, p);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STATIONARY DISTRIBUTION %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        %This will be the upperdiagonal of the matrix AAi
        AAupdiag        = zeros(p.Nb, 1); %This is necessary because of the peculiar way spdiags is defined.
        for j = 1:p.Na
            AAupdiag    =[AAupdiag; zeta(:, j, nz)];
        end

        %This will be the center diagonal of the matrix AAi
        AAcentdiag      = yy(:, 1, nz);
        for j = 2:p.Na-1
            AAcentdiag  = [AAcentdiag; yy(:, j, nz)];
        end
        AAcentdiag      = [AAcentdiag; yy(:, p.Na, nz)];

        %This will be the lower diagonal of the matrix AAi
        AAlowdiag       = chi(:, 2, nz);
        for j = 3:p.Na
            AAlowdiag   = [AAlowdiag; chi(:, j, nz)];
        end

        %Add up the upper, center, and lower diagonal into a sparse matrix
        AAi{nz}         = spdiags(AAcentdiag, 0, p.Nb * p.Na, p.Nb * p.Na) ...
            + spdiags(AAlowdiag, -p.Nb, p.Nb * p.Na, p.Nb * p.Na) ...
            + spdiags(AAupdiag, p.Nb, p.Nb * p.Na, p.Nb * p.Na);

    end

    AA                  = [AAi{1}, sparse(p.Nb * p.Na, p.Nb * p.Na); sparse(p.Nb * p.Na, p.Nb * p.Na), AAi{2}];
    A                   = AA + BB + Bswitch;

    M                   = p.Nb * p.Na * p.Nz;
    AT                  = A'; 
    % Fix one value so matrix isn't singular:
    vec                 = zeros(M, 1);
    iFix                = 1657;
    vec(iFix)           = .01;
    AT(iFix,:)          = [zeros(1, iFix - 1), 1, zeros(1, M - iFix)];

    % Solve system:
    g_tilde             = AT \ vec;
    g_sum               = g_tilde' * ones(p.Nb * p.Na * p.Nz, 1);
    g_tilde             = g_tilde ./ g_sum;
    
    g_stacked           = g_tilde ./ (p.dtildea_vec .* p.dtildeb_vec);
    g_stacked(g_stacked < 0)    = 0;

    g(:, :, 1)          = reshape(g_stacked(1:p.Nb * p.Na), p.Nb, p.Na);
    g(:, :, 2)          = reshape(g_stacked(p.Nb * p.Na + 1:p.Nb * p.Na * 2), p.Nb, p.Na);
    
    %% Pack
    sol.cpol            = c;
    sol.dpol            = d;
    sol.apol            = m;
    sol.bpol            = s;
    sol.sc              = sc;
    sol.sd              = sd;
    sol.dst             = g;
    sol.isvalid         = isvalid;
    sol.V               = V;
    
end