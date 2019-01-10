function d = two_asset_kinked_FOC(pa,pb,a,options,params)
    %% Unpack
    chi0    = params.chi0;
    chi1    = params.chi1;
    
    %% Return
    d       = min(pa./pb - 1 + chi0,0).*a/chi1 +  max(pa./pb - 1 - chi0,0).*a/chi1;
end