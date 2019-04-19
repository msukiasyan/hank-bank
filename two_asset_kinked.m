clear all; close all; clc;

%% Parameters
params.ga           = 2;                % CRRA utility with parameter gamma
params.rho          = 0.07;             % discount rate
params.chi0         = 0.03;              % 0.01;
params.chi1         = 2.0;
params.xi           = 0.01;              % fraction of income that is automatically deposited
params.Aprod        = 1.0;
params.alpha        = 0.4;
params.delta        = 0.05;
params.rho_bank     = -log(0.85);
params.f_bank       = 0.15;
params.theta_bank   = 0.4;

%Income process (two-state Poisson process):
params.Nz           = 2;
params.z            = [.4,1.7];          
params.la_mat       = [-1/5, 1/5; 1/5, -1/5];

params.crit         = 10^(-8);
params.Delta        = 2;
options.maxit       = 200;

params.r_minus      = 0.04;
params.r_F          = 0.042;
params              = get_all_params(options, params);

params.acurve       = 5;
params.bcurve       = 2;
params.Nb           = 60;
params.bmin         = 0;
params.bmax         = 40;
params.Na           = 100;
params.amin         = 0;
params.amax         = 60;

%% Options
options.debug_v     = 1;
options.debug_eq    = 0;

%% Setup
params              = setup(options, params);

%% Test
tic;
sol                 = get_policies(options, params);
toc
show_plots(options, params, sol, '-x');
return

%% Find equilibrium
rates               = fmincon(@(x) eqs(options, params, x), ...
    [0.045 0.045], [1, -1], 0, [], [], [0.0 0.0], [0.1 0.1], [], ...
    optimoptions('fmincon', 'Display', 'iter'));
return

%% Comparative statics
cs_theta_bank       = comp_stat(options, params, 'theta_bank', 0.38, 0.42, 5);

return

%% Plots

figure(5)
icut = 50; jcut=50;
bcut=params.b(1:icut); acut=params.a(1:jcut);
gcut = dst(1:icut,1:jcut,:);

set(gcf,'PaperPosition',[0 0 30 10])
subplot(2,2,3)
surf(bcut,acut,gcut(:,:,1)') %, 'LineStyle', 'none'
rotate3d on;
set(gca,'FontSize',16)
view([-10 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([params.bmin params.b(icut)])
ylim([params.amin params.a(jcut)])
title('Stationary, Low Type, Low Variance')

subplot(2,2,4)
surf(bcut,acut,gcut(:,:,2)') %, 'LineStyle', 'none'
rotate3d on;
set(gca,'FontSize',16)
view([-10 30])
xlabel('Liquid Wealth, b')
ylabel('Illiquid Wealth, a')
xlim([params.bmin params.b(icut)])
ylim([params.amin params.a(jcut)])
title('Stationary, High Type, Low Variance')


