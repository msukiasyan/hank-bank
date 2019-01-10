function eq = two_asset_kinked_cost(d,a,options,params)
    %% Unpack
    chi0    = params.chi0;
    chi1    = params.chi1;
    
    %% Return
    eq      = chi0.*abs(d) + chi1.*d.^2/2.*(max(a,10^(-5))).^(-1);
end