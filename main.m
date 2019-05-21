clear all; close all; clc;
warning off

%% Parameters
params.ga           = 2;                    % CRRA utility with parameter gamma
params.rho          = 0.01272892513;        % discount rate
params.chi0         = 0.04383;
params.chi1         = 0.95616994593;
params.chi2         = 1.40176;
params.chi3         = 0.03 *2.92/4.0;
params.xi           = 0.0000;               % fraction of income that is automatically deposited
params.Aprod        = 1.0;
params.alpha        = 0.33;
params.delta        = 0.07 / 4;
params.rho_bank     = -log(0.98);
params.f_bank       = 0.02;
params.theta_bank   = 0.3;

% Income process
params.Nz           = 3;
params.zfactor      = 1.0;
% params.zbase        = [0.8, 1.2] *  0.332349057201709;
% params.la_mat       = [-0.0667, 0.0667; 0.0667, -0.0667];
params.zbase        = [0.10, 0.50, 0.70];
params.zbase        = params.zbase / mean(params.zbase);
params.la_mat       = [-0.05, 0.05, 0.00; ...
                        0.03, -0.06, 0.03; ...
                        0.0, 0.05, -0.05];

params.dmax         = 1e10;
params.crit         = 10^(-8);
params.critKFE      = 10^(-12);
params.Delta        = 1000000.0;
params.DeltaKFE     = 1000.0;

params.r_minus      = 0.004;% 0.0380;
params.r_F          = 0.0045;%0.01422840843846188;% 0.0431;
params.r_init       = [params.r_minus, params.r_F];

%% Global parameters
glob                = struct();
params.acurve       = 1 / 0.15;
params.bcurve       = 1 / 0.35;
params.Nb           = 40;
params.bmin         = 0;
params.bmax         = 140;
params.Na           = 40;
params.amin         = 0;
params.amax         = 2000;
params.dtcurve      = 1 / 0.5;
params.Ndt          = 50;
params.dtmin        = 1 / 3;
params.dtmax        = 150;

%% Options
options.maxit       = 700;
options.maxitKFE    = 50000;
options.debug_v     = 1;
options.debug_eq    = 0;
options.stepK       = 0.1;
options.stepr       = 0.001;
options.maxittrans  = 500;
options.transtol    = 1e-6;

%% Setup
params              = setup(options, glob, params);

%% Test
% tic;
% params              = get_ss_params(options, glob, params);
% sol                 = get_policies(options, glob, params);
% toc
% stats               = calc_stats(options, glob, params, sol);
% 
% show_plots_ss(options, glob, params, stats);
% return

%% Asset supply
% Nra         = 5;
% Nrb         = 10;
% gridra      = linspace(0.0065, 0.0085, Nra);
% gridrb      = linspace(0.004, 0.006, Nrb);
% TSr         = zeros(Nra, Nrb);
% TDr         = zeros(Nra, Nrb);
% 
% parfor nr   = 1:(Nra * Nrb)
%     [nra, nrb]      = ind2sub([Nra, Nrb], nr);
%     par_copy        = params;
%     par_copy.r_F    = gridra(nra);
%     par_copy.r_plus = gridrb(nrb);
%     sol             = get_policies(options, par_copy);
%     stats           = calc_stats(options, par_copy, sol);
%     TSr(nr)         = stats.TS;
%     TDr(nr)         = stats.TD;
% end
% 
% figure;
% sgtitle('Asset Supply');
% 
% subplot(2, 2, 1);
% plot(gridrb, TSr);
% legend(string(1:Nra));
% xlabel('rb');
% ylabel('Illiquid Assets');
% title('TS(rb; ra)');
% 
% subplot(2, 2, 2);
% plot(gridra, TSr');
% legend(string(1:Nrb));
% xlabel('ra');
% ylabel('Illiquid Assets');
% title('TS(ra; rb)');
% 
% subplot(2, 2, 3);
% plot(gridrb, TDr);
% legend(string(1:Nra));
% xlabel('rb');
% ylabel('Liquid Assets');
% title('TD(rb; ra)');
% 
% subplot(2, 2, 4);
% plot(gridra, TDr');
% legend(string(1:Nrb));
% xlabel('ra');
% ylabel('Liquid Assets');
% title('TD(ra; rb)');
% return

%% Find equilibrium
% tic;
% [sol, stats]        = find_ss(options, glob, params, [], 1);
% toc;
% return

%% MIT shock to illiquid stock
tic;
[sol, stats]        = find_ss(options, glob, params, [], 1);
toc;
tic;
% paths               = transition_Nshock(options, glob, params, sol, stats, stats.NW * 0.95);
toc;
return

%% Comparative statics
tic;
cs_zfactor          = comp_stat(options, glob, params, 'zfactor', 0.95, 1.05, 5);
% cs_theta_bank       = comp_stat(options, glob, params, 'theta_bank', 0.48, 0.52, 5);
% cs_f_bank           = comp_stat(options, glob, params, 'f_bank', params.f_bank * 0.95, params.f_bank * 1.05, 5);
% cs_rho              = comp_stat(options, glob, params, 'rho', params.rho * 0.95, params.rho * 1.05, 5);
% cs_rho_bank         = comp_stat(options, glob, params, 'rho_bank', params.rho_bank * 0.95, params.rho_bank * 1.05, 5);
% cs_chi0             = comp_stat(options, glob, params, 'chi0', params.chi0 * 0.95, params.chi0 * 1.05, 5);
% cs_Aprod            = comp_stat(options, glob, params, 'Aprod', params.Aprod * 0.95, params.Aprod * 1.05, 5);
toc;

% load comp_stat.mat
% save comp_stat.mat
return

