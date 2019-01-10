function rs = eqs(options, params, inp)
    params.r_minus  = inp(1);
    params.r_F      = inp(2);
    params  = get_all_params(options, params);
    [cpol, dpol, bpol, apol, dst] = get_policies(options, params);

    ldist   = zeros(params.Nz, 1);
    ldist   = sum(dst(:, :, 1), 2) * params.da;

    ildist  = zeros(params.Nz, 1);
    ildist  = sum(dst(:, :, 1), 1) * params.db;

    TD      = ldist' * max(params.b, 0) * params.db;        % Total liquid deposits
    TB      = ldist' * max(-params.b, 0) * params.db;       % Total liquid borrowing
    TS      = ildist * params.a' * params.da;               % Total illiquid assets
    NW      = TS;                                           % Net worth = Total illiquid assets

    rs(1)   = TB - params.x_a * NW + params.K;
    rs(2)   = TD - params.x_a * NW + NW;
end