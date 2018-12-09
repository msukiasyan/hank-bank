function rs = eqs(pars)
    global Nz a b da db

    [cpol, dpol, bpol, apol, dst] = get_policies(pars(1), pars(2));
    TL = 0;
    TK = 0;

    for nz = 1:Nz
        TL = TL + sum(dst(:, :, nz), 2)' * b * db * da;
    end

    for nz = 1:Nz
        TK = TK + sum(dst(:, :, nz), 1) * a' * da * db;
    end

    rs(1) = TL;
    rs(2) = TK - pars(1);
end