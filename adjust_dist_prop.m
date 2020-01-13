function g1 = adjust_dist_prop(opt, glob, p, xgrid, xdelta, g, factor)
    gdist           = g / sum(g .* xdelta);
    g1cum           = interp1(xgrid, cumsum(gdist .* xdelta), xgrid / factor, 'linear', 'extrap');
    g1cum(end)      = 1.0;
    g1cum           = min(g1cum, 1.0);
    g1dist          = g1cum;
    g1dist(1)       = g1cum(1) / xdelta(1);
    g1dist(2:end)   = (g1cum(2:end) - g1cum(1:end-1)) ./ xdelta(2:end);
    g1dist(abs(g1dist) < 1e-15)     = 0.0;
    
    gmean           = sum(xgrid .* xdelta .* gdist);
    g1mean          = sum(xgrid .* xdelta .* g1dist);
    
    g1              = factor * gmean * g1dist / g1mean;
    if abs(g1mean)  < 1e-15
        g1          = g1dist * 0;
    end
    g1(1)           = (1.0 - sum(g1(2:end) .* xdelta(2:end))) / xdelta(1);
    
    g1              = g1 * sum(g .* xdelta);
end