function eq = adjustment_cost(d, a, opt, glob, p)
    %% Unpack
    chi0    = p.chi0;
    chi1    = p.chi1;
    chi2    = p.chi2;
    
    %% Return
    lincost = chi0 .* abs(d);
    eq      = lincost + chi1 .* abs(d ./ max(a,p.chi3)) .^ chi2 .* max(a,p.chi3);
end