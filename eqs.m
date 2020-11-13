function rs = eqs(opt, glob, p, inp)
    p.r_minus   = inp(1);
    p.r_F       = inp(2);
    p           = get_ss_params(opt, glob, p);
    if size(inp, 2) > 2
        p.r_plus    = inp(3);
    end
    p           = setup(opt, glob, p);
    if isfield(p, "NW") && p.NW < 0
        rs      = [1e10, 1e10];
        return
    end
%     p.r_bankeq  = p.r_F;
%     p.r_F       = p.mu * p.r_plus + (1 - p.mu) * p.r_F;
      sol         = get_policies(opt, glob, p);
    
    if ~sol.isvalid
        rs = 1e12;
        return
    end

    TD              = 0;
    TB              = 0;
    TS              = 0;
    p.dist_vec      = reshape(sol.dst, p.Nb * p.Na * p.Nz, 1);
    wt              = p.dtildea_vec .* p.dtildeb_vec .* p.dist_vec;
    h_vec           = reshape(sol.hpol, p.Nb * p.Na * p.Nz, 1);
    z_vec           = reshape(p.zzz, p.Nb * p.Na * p.Nz, 1);
    p.H             = sum(wt .* z_vec .* h_vec);   % labor supply
    p.K             = p.H * p.K_H;        % get capital
    
    ldist           = cell(p.Nz,1); 
    ildist          = cell(p.Nz,1);
    for nz = 1:p.Nz
        ldist{nz}   = trapz(p.a, sol.dst(:, :, nz), 2);         % sum(sol.dst(:, :, nz), 2) * p.da;
        ildist{nz}  = trapz(p.b, sol.dst(:, :, nz), 1);         % sum(sol.dst(:, :, nz), 1) * p.db;
        
        TD          = TD + trapz(max(p.b, 0), ldist{nz} .* p.b);       % Total liquid deposits
        TB          = TB + trapz(max(-p.b, 0), ldist{nz} .* p.b);      % Total liquid borrowing
        TS          = TS + trapz(p.a, ildist{nz} .* p.a);              % Total illiquid assets
    end

    NW              = TS;                                              % Net worth = Illiquid
    TD_bank         = (p.x_a - 1) * NW ;

    rt(1)           = TB - p.x_a * NW + p.K ;
    rt(2)           = TD - TD_bank;
    rs              = rt;
    

    if opt.debug_eq
        disp(rt(1));
        disp(rt(2));
    end
end