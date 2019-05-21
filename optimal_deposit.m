function d = optimal_deposit(pa, pb, a, opt,  glob, p)
    d       = (-(-min(pa./pb - 1 + p.chi0,0)) .^ (1 / (p.chi2 - 1)) .* max(a,p.chi3) / ((p.chi1 * p.chi2) .^ (1 / (p.chi2 - 1))) + ...
                max(pa./pb - 1 - p.chi0,0) .^ (1 / (p.chi2 - 1)) .* max(a,p.chi3) / ((p.chi1 * p.chi2) .^ (1 / (p.chi2 - 1))));
    d       = min(d, p.dmax);
end