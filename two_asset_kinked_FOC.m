function d = two_asset_kinked_FOC(pa,pb,a,options,params)
    %% Unpack
    chi0    = params.chi0;
    chi1    = params.chi1;
    chi2    = params.chi2;
    
    %% Return
    d       = zeros(size(pa));
    slct    = pa./pb - 1 + chi0 < 0;
    val     = abs(pa./pb - 1 + chi0) .^ (1 / (chi2 - 1)) / ((chi1 * chi2) .^ (1 / (chi2 - 1)));
    d(slct)     = -val(slct);
    d(~slct)    = val(~slct);
    
    d       = d .* max(a,params.chi3);
    d       = min(d, params.dmax);
%     d       = (-(-min(pa./pb - 1 + chi0,0)) .^ (1 / (chi2 - 1)) .* max(a,params.chi3) / ((chi1 * chi2) .^ (1 / (chi2 - 1))) + ...
%                 max(pa./pb - 1 - chi0,0) .^ (1 / (chi2 - 1)) .* max(a,params.chi3) / ((chi1 * chi2) .^ (1 / (chi2 - 1))));
%     d       = min(d, params.dmax);
end