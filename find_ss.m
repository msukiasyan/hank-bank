function [sol, stats] = find_ss(opt, glob, p, iguess, pl)
    if nargin < 3
        iguess  = p.r_init;
        pl      = 0;
    elseif nargin < 4
        pl      = 0;
    end
    if isempty(iguess)
        iguess  = p.r_init;
    end
    
    % rates               = fmincon(@(x) eqs(options, params, x), ...
    %     [0.038 0.043], [1, -1], 0, [], [], [0.0 0.0], [0.1 0.1], [], ...
    %     optimoptions('fmincon', 'Display', 'iter'));
    p                   = get_all_params(opt, p);
    rates               = fsolve(@(x) eqs(opt, p, x), ...
        iguess, optimoptions('fsolve', 'Display', 'iter','UseParallel', true));
    p.r_F               = rates(2);
    p.r_minus           = rates(1);
    p                   = get_all_params(opt, p);
    p                   = setup(opt, p);
    sol                 = get_policies(opt, p);
    stats               = calc_stats(opt, p, sol);
    if pl == 1
        show_plots_ss(opt, p, stats);
    end
end