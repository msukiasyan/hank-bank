clear all; close all; clc;
warning off

%% Parameters
params.ga           = 2;                % CRRA utility with parameter gamma
params.rho          = 0.01272892513;             % discount rate
params.chi0         = 0.04383;              % 0.01;
params.chi1         = 0.95616994593;
params.chi2         = 1.40176;
params.chi3         = 0.03 *2.92/4.0;
params.xi           = 0.0000;              % fraction of income that is automatically deposited
params.Aprod        = 1.0;
params.alpha        = 0.33;
params.delta        = 0.07 / 4;
params.rho_bank     = -log(0.98);
params.f_bank       = 0.025;
params.theta_bank   = 0.5;

%Income process (two-state Poisson process):
params.Nz           = 3;
params.zfactor      = 1.0;
% params.zbase        = [0.8, 1.2] *  0.332349057201709;
% params.la_mat       = [-0.0667, 0.0667; 0.0667, -0.0667];
params.zbase        = [0.10, 0.50, 0.70];
params.la_mat       = [-0.05, 0.05, 0.00; ...
                        0.03, -0.06, 0.03; ...
                        0.0, 0.05, -0.05];
params.la_mat       = params.la_mat * 1.00;

params.dmax         = 1e10;
params.crit         = 10^(-8);
params.critKFE      = 10^(-12);
params.Delta        = 1000000.0;
params.DeltaKFE     = 1000.0;
options.maxit       = 700;
options.maxitKFE    = 50000;

params.r_minus      = 0.0070;% 0.0380;
params.r_F          = 0.0082;%0.01422840843846188;% 0.0431;
params.r_init       = [params.r_minus, params.r_F];
params              = get_all_params(options, params);

%% %%%%%%%%%%%%%%%%%
params.r_plus       = 0.005;
params.w            = 1.27005019392697;
%% %%%%%%%%%%%%%%%%%

params.acurve       = 1 / 0.15; % 1.5;
params.bcurve       = 1 / 0.35; % 2.0;
params.Nb           = 40;
params.bmin         = 0;
params.bmax         = 140;
params.Na           = 40;
params.amin         = 0;
params.amax         = 2000;

%% Options
options.debug_v     = 1;
options.debug_eq    = 0;

%% Setup
params              = setup(options, params);

%% Test
tic;
sol                 = get_policies(options, params);
toc
stats               = calc_stats(options, params, sol);
show_plots_ss(options, params, stats);
return

%% Find equilibrium
tic;
[sol, stats]        = find_ss(options, params, [], 1);
toc;
return

%% Comparative statics
tic;
cs_zfactor          = comp_stat(options, params, 'zfactor', 0.95, 1.05, 5);
cs_theta_bank       = comp_stat(options, params, 'theta_bank', 0.38, 0.42, 5);
cs_f_bank           = comp_stat(options, params, 'f_bank', 0.13, 0.17, 5);
toc;

save comp_stat.mat
return

