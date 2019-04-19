function show_plots(opt, p, s, stl)
    if nargin < 4
        stl = '-';
    end
    TD              = 0;
    TB              = 0;
    TS              = 0;
    ldist           = cell(p.Nz,1);
    ildist          = cell(p.Nz,1);
    for nz = 1:p.Nz
        ldist{nz}   = trapz(p.a, s.dst(:, :, nz), 2);         
        ildist{nz}  = trapz(p.b, s.dst(:, :, nz), 1);         
        
        TD          = TD + trapz(max(p.b, 0), ldist{nz} .* p.b);       % Total liquid deposits
        TB          = TB + trapz(max(-p.b, 0), ldist{nz} .* p.b);      % Total liquid borrowing
        TS          = TS + trapz(p.a, ildist{nz} .* p.a);              % Total illiquid assets
    end

    NW              = TS;                                       % Net worth = Total illiquid assets
    
    figure;
    subplot(2,2,1);
    plot(p.b, ldist{1}, stl)

    subplot(2,2,2);
    plot(p.a, ildist{1}, stl)

    subplot(2,2,3);
    plot(p.b, ldist{2}, stl)

    subplot(2,2,4);
    plot(p.a, ildist{2}, stl)
end