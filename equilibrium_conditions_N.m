% Unpack vars

V =
alpha_hat = 
g = 
N =

logAggregateTFP =



K =
r_plus =
alpha_hat = 
lambda =
q_K =

w = exp(logAggregateTFP) * (1 - alpha) * (K ^ alpha); 
r_minus = exp(logAggregateTFP) * alpha * (K ^ (alpha -1) / q_K ) + q_K_dot/q_K -delta; 

leverage = alpha_hat/theta;

lambda = leverage*(r_minus-r_plus);
residual_alpha_hat = (rho_B + f)*alpha_hat - f - alpha_hat*(r_plus-lambda) - alpha_hat_dot;




consumptionEError = vars(2*nVars+1,1);
inflationEError = vars(2*nVars+2,1);

rhoSHOCK = vars(2*nVars+3);
mp_shockSHOCK = vars(2*nVars+4);

r_Nominal = rho_bar + phi * inflation + mp_shock;
r = r_Nominal - inflation;

L_y = (Y/Z_y)^(1/alpha_y);
L_n= (N/Z_n)^(1/alpha_n);







%% HJB

%% Collect dynamics
% HJB equation
EulerResidual = consumptionDot + consumptionEError - consumption*(r-rho)/gamma;

% Inflation
pcResidual = (rho - consumptionDot/consumption) * inflation + varieties^(-1) *(sigma-1)/kappa * (1-sigma/(sigma-1)/mu) - inflationDot + inflationEError;
%pcResidual = -((r - 0) * inflation - (sigma/kappa * (Pw_P- (sigma-1)/sigma) + inflationDot - inflationEError));

% Exog. variables
MPResidual = mp_shock_Dot - (-mp_shock_lag * mp_shock + ssigma_mp_shock * mp_shockSHOCK);
rhoResidual = rho_Dot - (-rho_lag * (rho-rho_bar) + ssigma_rho * rhoSHOCK);

% Market Clearing
FocLnResidual = N*L_y/L_n - 1/(mu-1) * alpha_y/alpha_n;
marketResidual = alpha_y*Y/L_y - (consumption)^gamma * mu * disutil * (L_n + N*L_y)^varphi;
consresidual = consumption - N*Y;


v_residual = [EulerResidual;
              pcResidual;
              rhoResidual;
              MPResidual;
              FocLnResidual;
              marketResidual;
              consresidual;
];
