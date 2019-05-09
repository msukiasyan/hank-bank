function eq = two_asset_kinked_cost(d,a,options,params)
    %% Unpack
    chi0    = params.chi0;
    chi1    = params.chi1;
    chi2    = params.chi2;
    
    %% Return
    if d > 0
        lincost = chi0.*abs(d);
    else
        lincost = chi0.*abs(d);
    end
    eq      = lincost + chi1 .* abs(d ./ max(a,params.chi3)) .^ chi2 .* max(a,params.chi3);
end