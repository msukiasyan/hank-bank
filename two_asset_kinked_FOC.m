function d = two_asset_kinked_FOC(pa,pb,a)
global chi0 chi1
d = min(pa./pb - 1 + chi0,0).*a/chi1 +  max(pa./pb - 1 - chi0,0).*a/chi1;