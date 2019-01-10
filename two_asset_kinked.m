clear all; close all; clc;

%% Setup
params.ga           = 2;                % CRRA utility with parameter gamma
params.rho          = 0.07;             % discount rate
params.chi0         = 0.1;              % 0.01;
params.chi1         = 2.0;
params.xi           = 0.1;              % fraction of income that is automatically deposited
params.Aprod        = 1.0;
params.alp          = 0.4;
params.delt         = 0.05;
params.rho_bank     = -log(0.96);
params.f_bank       = 0.05;
params.theta_bank   = 0.19;

%Income process (two-state Poisson process):
params.Nz           = 2;
params.z            = [.5,1.5];          
params.la_mat       = [-1/2, 1/2; 1/2, -1/2];

options.crit        = 10^(-5);
options.Delta       = 100;
options.maxit       = 35;

params.r_minus      = 0.04;
params.r_F          = 0.05;
params              = get_all_params(options, params);

% Grids
params.I            = 101;
params.bmin         = -2;
params.bmax         = 48;
params.b            = linspace(params.bmin,params.bmax,params.I)';
params.db           = (params.bmax-params.bmin)/(params.I-1);

params.J            = 100;
params.amin         = 0;
params.amax         = 150;
params.a            = linspace(params.amin,params.amax,params.J);
params.da           = (params.amax-params.amin)/(params.J-1);

%% Test

% [cpol, dpol, bpol, apol, dst] = get_policies(0.05, 0.04, 0.0382, options, params);
% 
% ldist   = zeros(params.Nz, 1);
% ldist   = sum(dst(:, :, 1), 2) * params.da;
% 
% ildist  = zeros(params.Nz, 1);
% ildist  = sum(dst(:, :, 1), 1) * params.db;
% 
% TD      = ldist' * max(params.b, 0) * params.db;        % Total liquid deposits
% TB      = ldist' * max(-params.b, 0) * params.db;       % Total liquid borrowing
% TS      = ildist * params.a' * params.da;               % Total illiquid assets
% NW      = TS;                                           % Net worth = Total illiquid assets
% 
% 
% subplot(1,2,1);
% plot(params.b, ldist)
% 
% subplot(1,2,2);
% plot(params.a, ildist)
% return

sol = fsolve(@(x) eqs(options, params, x), [0.04, 0.06], optimoptions('fsolve', 'Display', 'iter'))
return

%% Plots

figure(5)
icut = 50; jcut=50;
bcut=b(1:icut); acut=a(1:jcut);
gcut = dst(1:icut,1:jcut,:);

set(gcf,'PaperPosition',[0 0 30 10])
subplot(2,2,3)
surf(bcut,acut,gcut(:,:,1)') %, 'LineStyle', 'none'
rotate3d on;
set(gca,'FontSize',16)
view([-10 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin b(icut)])
ylim([amin a(jcut)])
title('Stationary, Low Type, Low Variance')

subplot(2,2,4)
surf(bcut,acut,gcut(:,:,2)') %, 'LineStyle', 'none'
rotate3d on;
set(gca,'FontSize',16)
view([-10 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([bmin b(icut)])
ylim([amin a(jcut)])
title('Stationary, High Type, Low Variance')


