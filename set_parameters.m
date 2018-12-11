global alpha delta xi a1 theta f rho rhoB M_N ga
alpha = 0.35; % Production function F = K^alpha * L^(1-alpha) 
delta = 0.1;  % Capital depreciation
xi = 0.5; % adjustment cost elasticity
a1 = 1/(1-xi)*delta^xi; % adjustment cost parameter 1
theta = 0.5; %fraction of divertable assets 
f = 0.1; % exit rate of banks
rho = 0.05;   % discount rate HH
rhoB = 0.05;   % discount rate BANK
M_N = 0.01; % recapitalization rate
ga = 1; % leave 1 for log utility