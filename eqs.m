function rs = eqs(opt, p, inp, pl)
    if nargin < 4
        pl      = 0;
    end

    p.r_minus   = inp(1);
    p.r_F       = inp(2);
    p           = get_all_params(opt, p);
    if size(inp, 2) > 2
        p.r_plus    = inp(3);
    end
    sol         = get_policies(opt, p);
    
    if ~sol.isvalid
        rs = 1e12;
        return
    end

    TD              = 0;
    TB              = 0;
    TS              = 0;
    ldist           = cell(p.Nz,1);
    ildist          = cell(p.Nz,1);
    for nz = 1:p.Nz
        ldist{nz}   = trapz(p.a, sol.dst(:, :, nz), 2);         % sum(sol.dst(:, :, nz), 2) * p.da;
        ildist{nz}  = trapz(p.b, sol.dst(:, :, nz), 1);         % sum(sol.dst(:, :, nz), 1) * p.db;
        
        TD          = TD + trapz(max(p.b, 0), ldist{nz} .* p.b);       % Total liquid deposits
        TB          = TB + trapz(max(-p.b, 0), ldist{nz} .* p.b);      % Total liquid borrowing
        TS          = TS + trapz(p.a, ildist{nz} .* p.a);              % Total illiquid assets
    end

    NW              = TS;                                       % Net worth = Total illiquid assets
    
    if pl == 1
        show_plots(opt, p, sol, '-x');
    end

    rt(1)           = TB - p.x_a * NW + p.K;
    rt(2)           = TD - p.x_a * NW + NW;
    rs              = rt(1) ^ 2 + rt(2) ^ 2;
    
    if opt.debug_eq
        disp(rt(1));
        disp(rt(2));
    end
end