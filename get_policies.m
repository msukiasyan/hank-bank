function [cpol, dpol, bpol, apol, dst] = get_policies(options, params)

    %% Unpack
    chi0        = params.chi0;
    chi1        = params.chi1;
    ga          = params.ga;
    rho         = params.rho; 
    xi          = params.xi;
    Aprod       = params.Aprod;
    alp         = params.alp;
    delt        = params.delt;
    I           = params.I;
    J           = params.J;
    da          = params.da;
    db          = params.db;
    a           = params.a;
    b           = params.b;
    Nz          = params.Nz;
    bmin        = params.bmin;
    amax        = params.amax;
    z           = params.z;
    la_mat      = params.la_mat;
    w           = params.w;
    bmax        = params.bmax;
    amin        = params.amin;
    
    crit        = options.crit;
    Delta       = options.Delta;
    maxit       = options.maxit;
    
    %% Solve
    
    ra          = params.r_F;
    rb_pos      = params.r_plus;
    rb_neg      = params.r_minus;

    if ra - 1/chi1 > 0
        disp('Warning: ra - 1/chi1 > 0')
    end

    bb = b*ones(1,J);
    aa = ones(I,1)*a;
    zz = ones(J,1)*z;

    bbb = zeros(I,J,Nz); aaa = zeros(I,J,Nz); zzz = zeros(I,J,Nz);
    for nz=1:Nz
        bbb(:,:,nz) = bb;
        aaa(:,:,nz) = aa;
        zzz(:,:,nz) = z(nz);
    end

    Bswitch = [
        speye(I*J)*la_mat(1,1), speye(I*J)*la_mat(1,2);
        speye(I*J)*la_mat(2,1), speye(I*J)*la_mat(2,2)];

    %Preallocation
    VbF = zeros(I,J,Nz);
    VbB = zeros(I,J,Nz);
    VaF = zeros(I,J,Nz);
    VaB = zeros(I,J,Nz);
    c = zeros(I,J,Nz);
    updiag = zeros(I*J,Nz);
    lowdiag = zeros(I*J,Nz);
    centdiag = zeros(I*J,Nz);
    AAi = cell(Nz,1);
    BBi = cell(Nz,1);


    %INITIAL GUESS
    v0 = (((1-xi)*w*zzz + ra.*aaa + 0.06 .*bbb).^(1-ga))/(1-ga)/rho;
    v = v0;

    %return at different points in state space
    %matrix of liquid returns
    Rb = rb_pos.*(bbb>0) + rb_neg.*(bbb<0);
    raa = ra.*ones(1,J);
    %if ra>>rb, impose tax on ra*a at high a, otherwise some households
    %accumulate infinite illiquid wealth (not needed if ra is close to or less than rb)
    tau = 10; raa = ra.*(1 - (1.33.*amax./a).^(1-tau));
    %matrix of illiquid returns
    Ra(:,:,1) = ones(I,1)*raa;
    Ra(:,:,2) = ones(I,1)*raa;


    for n=1:maxit
        V = v;   
        %DERIVATIVES W.R.T. b
        % forward difference
        VbF(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
        VbF(I,:,:) = ((1-xi)*w*zzz(I,:,:) + Rb(I,:,:).*bmax).^(-ga); %state constraint boundary condition
        % backward difference
        VbB(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/db;
        VbB(1,:,:) = ((1-xi)*w*zzz(1,:,:) + Rb(1,:,:).*bmin).^(-ga); %state constraint boundary condition

        %DERIVATIVES W.R.T. a
        % forward difference
        VaF(:,1:J-1,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da;
        % backward difference
        VaB(:,2:J,:) = (V(:,2:J,:)-V(:,1:J-1,:))/da;

        %useful quantities
        c_B = max(VbB,10^(-6)).^(-1/ga);
        c_F = max(VbF,10^(-6)).^(-1/ga);
        dBB = two_asset_kinked_FOC(VaB,VbB,aaa,options,params);
        dFB = two_asset_kinked_FOC(VaB,VbF,aaa,options,params);
        %VaF(:,J,:) = VbB(:,J,:).*(1-ra.*chi1 - chi1*w*zzz(:,J,:)./a(:,J,:));
        dBF = two_asset_kinked_FOC(VaF,VbB,aaa,options,params);
        %VaF(:,J,:) = VbF(:,J,:).*(1-ra.*chi1 - chi1*w*zzz(:,J,:)./a(:,J,:));
        dFF = two_asset_kinked_FOC(VaF,VbF,aaa,options,params);

        %UPWIND SCHEME
        d_B = (dBF>0).*dBF + (dBB<0).*dBB;
        %state constraints at amin and amax
        d_B(:,1,:) = (dBF(:,1,:)>10^(-12)).*dBF(:,1,:); %make sure d>=0 at amax, don't use VaB(:,1,:)
        d_B(:,J,:) = (dBB(:,J,:)<-10^(-12)).*dBB(:,J,:); %make sure d<=0 at amax, don't use VaF(:,J,:)
        d_B(1,1,:)=max(d_B(1,1,:),0);
        %split drift of b and upwind separately
        sc_B = (1-xi)*w*zzz + Rb.*bbb - c_B;
        sd_B = (-d_B - two_asset_kinked_cost(d_B,aaa,options,params));

        d_F = (dFF>0).*dFF + (dFB<0).*dFB;
        %state constraints at amin and amax
        d_F(:,1,:) = (dFF(:,1,:)>10^(-12)).*dFF(:,1,:); %make sure d>=0 at amin, don't use VaB(:,1,:)
        d_F(:,J,:) = (dFB(:,J,:)<-10^(-12)).*dFB(:,J,:); %make sure d<=0 at amax, don't use VaF(:,J,:)

        %split drift of b and upwind separately
        sc_F = (1-xi)*w*zzz + Rb.*bbb - c_F;
        sd_F = (-d_F - two_asset_kinked_cost(d_F,aaa,options,params));
        sd_F(I,:,:) = min(sd_F(I,:,:),0);

        Ic_B = (sc_B < -10^(-12));
        Ic_F = (sc_F > 10^(-12)).*(1- Ic_B);
        Ic_0 = 1 - Ic_F - Ic_B;

        Id_F = (sd_F > 10^(-12));
        Id_B = (sd_B < -10^(-12)).*(1- Id_F);
        Id_B(1,:,:)=0;
        Id_F(I,:,:) = 0; Id_B(I,:,:) = 1; %don't use VbF at bmax so as not to pick up articial state constraint
        Id_0 = 1 - Id_F - Id_B;

        c_0 = (1-xi)*w*zzz + Rb.*bbb;

        c = c_F.*Ic_F + c_B.*Ic_B + c_0.*Ic_0;
        u = c.^(1-ga)/(1-ga);

        %CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
        X = -Ic_B.*sc_B/db -Id_B.*sd_B/db;
        Y = (Ic_B.*sc_B - Ic_F.*sc_F)/db + (Id_B.*sd_B - Id_F.*sd_F)/db;
        Z = Ic_F.*sc_F/db + Id_F.*sd_F/db;

        for i = 1:Nz
            centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
        end

        lowdiag(1:I-1,:) = X(2:I,1,:);
        updiag(2:I,:) = Z(1:I-1,1,:);
        for j = 2:J
            lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
            updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
        end

        for nz=1:Nz
            BBi{nz}=spdiags(centdiag(:,nz),0,I*J,I*J)+spdiags([updiag(:,nz);0],1,I*J,I*J)+spdiags([lowdiag(:,nz);0],-1,I*J,I*J);
        end

        BB = [BBi{1}, sparse(I*J,I*J); sparse(I*J,I*J), BBi{2}];

        %CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
        dB = Id_B.*dBB + Id_F.*dFB;
        dF = Id_B.*dBF + Id_F.*dFF;
        MB = min(dB,0);
        MF = max(dF,0) + xi*w*zzz + Ra.*aaa;
        MB(:,J,:) = xi*w*zzz(:,J,:) + dB(:,J,:) + Ra(:,J,:).*amax; %this is hopefully negative
        MF(:,J,:) = 0;
        chi = -MB/da;
        yy =  (MB - MF)/da;
        zeta = MF/da;

        %MATRIX AAi
        for nz=1:Nz
            %This will be the upperdiagonal of the matrix AAi
            AAupdiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
            for j=1:J
                AAupdiag=[AAupdiag;zeta(:,j,nz)];
            end

            %This will be the center diagonal of the matrix AAi
            AAcentdiag= yy(:,1,nz);
            for j=2:J-1
                AAcentdiag=[AAcentdiag;yy(:,j,nz)];
            end
            AAcentdiag=[AAcentdiag;yy(:,J,nz)];

            %This will be the lower diagonal of the matrix AAi
            AAlowdiag=chi(:,2,nz);
            for j=3:J
                AAlowdiag=[AAlowdiag;chi(:,j,nz)];
            end

            %Add up the upper, center, and lower diagonal into a sparse matrix
            AAi{nz} = spdiags(AAcentdiag,0,I*J,I*J)+spdiags(AAlowdiag,-I,I*J,I*J)+spdiags(AAupdiag,I,I*J,I*J);

        end

        AA = [AAi{1}, sparse(I*J,I*J); sparse(I*J,I*J), AAi{2}];

        A = AA + BB + Bswitch;

        if max(abs(sum(A,2)))>10^(-9)
           disp('Improper Transition Matrix')
           break
        end

        B = (1/Delta + rho)*speye(I*J*Nz) - A;

        u_stacked = reshape(u,I*J*Nz,1);
        V_stacked = reshape(V,I*J*Nz,1);

        vec = u_stacked + V_stacked/Delta;

        V_stacked = B\vec; %SOLVE SYSTEM OF EQUATIONS

        V = reshape(V_stacked,I,J,Nz);   


        Vchange = V - v;
        v = V;

        b_dist(:,n) = max(abs(Vchange(:,:,2)),[],2);
        a_dist(:,n) = max(max(abs(Vchange),[],3),[],1);
        ab_dist(:,:,n) = max(abs(Vchange),[],3);

        dist(n) = max(max(max(abs(Vchange))));
        disp(['Value Function, Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
        if dist(n)<crit
            disp('Value Function Converged, Iteration = ')
            disp(n)
            break
        end

    end

    d = Id_B.*d_B + Id_F.*d_F;
    m = d + xi*w*zzz + Ra.*aaa;
    s = (1-xi)*w*zzz + Rb.*bbb - d - two_asset_kinked_cost(d,aaa,options,params) - c;

    sc = (1-xi)*w*zzz + Rb.*bbb - c;
    sd = - d - two_asset_kinked_cost(d,aaa,options,params);
    
    cpol = c;
    dpol = d;
    apol = m;
    bpol = s;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STATIONARY DISTRIBUTION %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RECONSTRUCT TRANSITION MATRIX WITH SIMPLER UPWIND SCHEME
    X = -min(s,0)/db;
    Y = min(s,0)/db - max(s,0)/db;
    Z = max(s,0)/db;

    for i = 1:Nz
        centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
    end

    lowdiag(1:I-1,:) = X(2:I,1,:);
    updiag(2:I,:) = Z(1:I-1,1,:);
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];
        updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end

    for nz=1:Nz
        BBi{nz}=spdiags(centdiag(:,nz),0,I*J,I*J)+spdiags([updiag(:,nz);0],1,I*J,I*J)+spdiags([lowdiag(:,nz);0],-1,I*J,I*J);
    end

    BB = [BBi{1}, sparse(I*J,I*J); sparse(I*J,I*J), BBi{2}];

    %CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
    chi = -min(m,0)/da;
    yy =  min(m,0)/da - max(m,0)/da;
    zeta = max(m,0)/da;


    %MATRIX AAi
    for nz=1:Nz
        %This will be the upperdiagonal of the matrix AAi
        AAupdiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
        for j=1:J
            AAupdiag=[AAupdiag;zeta(:,j,nz)];
        end

        %This will be the center diagonal of the matrix AAi
        AAcentdiag= yy(:,1,nz);
        for j=2:J-1
            AAcentdiag=[AAcentdiag;yy(:,j,nz)];
        end
        AAcentdiag=[AAcentdiag;yy(:,J,nz)];

        %This will be the lower diagonal of the matrix AAi
        AAlowdiag=chi(:,2,nz);
        for j=3:J
            AAlowdiag=[AAlowdiag;chi(:,j,nz)];
        end

        %Add up the upper, center, and lower diagonal into a sparse matrix
        AAi{nz} = spdiags(AAcentdiag,0,I*J,I*J)+spdiags(AAlowdiag,-I,I*J,I*J)+spdiags(AAupdiag,I,I*J,I*J);

    end

    AA = [AAi{1}, sparse(I*J,I*J); sparse(I*J,I*J), AAi{2}];
    A = AA + BB + Bswitch;

    M = I*J*Nz;
    AT = A';
    % Fix one value so matrix isn't singular:
    vec = zeros(M,1);
    iFix = 1657;
    vec(iFix)=.01;
    AT(iFix,:) = [zeros(1,iFix-1),1,zeros(1,M-iFix)];

    % Solve system:
    g_stacked = AT\vec;
    g_sum = g_stacked'*ones(M,1)*da*db;
    g_stacked = g_stacked./g_sum;


    g(:,:,1) = reshape(g_stacked(1:I*J),I,J);
    g(:,:,2) = reshape(g_stacked(I*J+1:I*J*2),I,J);
    dst = g;
    
end