function res = prctilew(x, w, p)
    [x_sorted, x_order]         = sort(x);
    w_sorted                    = w(x_order);
    w_cum                       = cumsum(w_sorted) - 0.5 * w_sorted;
    w_cum                       = w_cum / w_cum(end);
    w_cum                       = cumsum(ones(size(w_cum))) * 2 * eps + w_cum;
    [w_unique, ind]             = unique(w_cum);
    res                         = interp1(w_unique, x_sorted(ind), p, 'linear', 'extrap');
end