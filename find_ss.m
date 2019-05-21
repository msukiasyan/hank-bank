function [sol, stats] = find_ss(opt, glob, p, iguess, pl)
    if nargin < 4
        iguess  = p.r_init;
        pl      = 0;
    elseif nargin < 5
        pl      = 0;
    end
    if isempty(iguess)
        iguess  = p.r_init;
    end
    
    p                   = get_ss_params(opt, glob, p);
    rates               = fsolve(@(x) eqs(opt, glob, p, x), ...
        iguess, optimoptions('fsolve', 'Display', 'iter','UseParallel', true));
    p.r_F               = rates(2);
    p.r_minus           = rates(1);
    p                   = get_ss_params(opt, glob, p);
    sol                 = get_policies(opt, glob, p);
    stats               = calc_stats(opt, glob, p, sol);
    if pl == 1
        show_plots_ss(opt, glob, p, stats);
    end
end