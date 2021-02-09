clear all; clc;
warning off

%% Parameters
params.ga           = 2;                    % CRRA utility with parameter gamma
params.rho          = 0.01272892513;        % discount rate
% params.chi0         = 1e5;
params.chi0         = 0.04383;
params.chi1         = 0.95616994593;
params.chi2         = 1.40176;
params.chi3         = 0.03 *2.92/4.0;
params.xi           = 0.0000;               % fraction of income that is automatically deposited
params.Aprod        = 1.0;
params.alpha        = 0.33;
params.delta        = 0.07 / 4;
params.rho_bank     = 0.02;
params.f_bank       = 0.04;
params.theta_bank   = 0.6;
params.kappa        = 1;
params.mu           = 0.0;                  % fraction of illiquid assets held in "illiquid deposits"
params.mu_bank      = 0.0;                  % fraction of bank deposits being illiquid
params.disutil      = 6.0;                  % labor disutility parameter
params.frisch       = 1.0;                  % Frisch elasticity of labor supply


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
params.crit         = 10^(-9);
params.critKFE      = 10^(-12);
params.Delta        = 1000000.0;
params.DeltaKFE     = 1000.0;

params.r_minus      = 0.004;
params.r_F          = 0.0062;%0.01422840843846188;
params.r_init       = [params.r_minus, params.r_F];

%% Global parameters
glob                = struct();
params.acurve       = 1 / 0.15;
params.bcurve       = 1 / 0.35;
params.Nb           = 40;
params.bmin         = 0;
params.bmax         = 150;
params.Na           = 40;
params.amin         = 0;
params.amax         = 2000;
params.dtcurve      = 1 / 0.5;
params.Ndt          = 40;
params.dtmin        = 1 / 5;
params.dtmax        = 100;

%% GK (family) parameters
params.distGK       = "twopoint";
params.fracGK       = 1.00;          % Fraction of people having bank ownership

params.MGK          = 0.1;

%% Options
options.maxit       = 2000;
options.maxitKFE    = 50000;
options.debug_v     = 1;
options.debug_eq    = 1;
options.stepK       = 0.1;
options.stepr       = 0.001;
options.maxittrans  = 1000;
options.transtol    = 1e-6;
options.divtoliq    = true;

options.GK          = 0;                         % If 1, must also be that divtoliq == true and chi0 large

options.stepK_nt    = 1e-4;
options.stepx_a_nt  = 1e-5;
options.stepH_nt    = 1e-5;
options.tol_nt      = 1e-10;
options.maxittrans_nt   = 140;
options.debug_trans = 1;

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

%% Find steady state equilibrium
% tic;
% [sol, stats]        = find_ss(options, glob, params, [], 1);
% toc;
% return

%% MIT shock to illiquid stock
% tic;
% [sol, stats]        = find_ss(options, glob, params, [], 0);
% toc;
% tic;
% [paths, statst]     = transition_Nshock_newton(options, glob, params, sol, stats, stats.NW * 0.95);
% show_plots_mit(options, glob, params, stats, paths, statst);
% toc;
% return

%% MIT shock
tic;
[sol, stats]        = find_ss(options, glob, params, [], 0);
toc;
tic;
[paths, statst]     = transition_Ashock_newton(options, glob, params, sol, stats, -0.01, 0.05,0);
[paths_benchmark, statst_benchmark]     = transition_Ashock_newton(options, glob, params, sol, stats, -0.01, 0.05,1);
[statst]            = calc_distr_irf(options, glob, params, statst,statst_benchmark);

toc;

show_plots_mit(options, glob, params, stats, paths, statst);
show_plots_distributional(options, glob, params, stats, paths, statst);

return

%% Comparative statics
% tic;
% % cs_fracGK           = comp_stat(options, glob, params, 'fracGK', 0.05, 0.95, 10);
% % cs_theta_bank       = comp_stat(options, glob, params, 'theta_bank', 0.25, 0.35, 5);
% % cs_zfactor          = comp_stat(options, glob, params, 'zfactor', 0.9, 1.1, 5);
% % cs_f_bank           = comp_stat(options, glob, params, 'f_bank', params.f_bank * 0.95, params.f_bank * 1.05, 5);
% % cs_rho              = comp_stat(options, glob, params, 'rho', params.rho * 0.95, params.rho * 1.05, 5);
% % cs_rho_bank         = comp_stat(options, glob, params, 'rho_bank', params.rho_bank * 0.95, params.rho_bank * 1.05, 5);
% % cs_chi0             = comp_stat(options, glob, params, 'chi0', params.chi0 * 0.95, params.chi0 * 1.05, 5);
% % cs_Aprod            = comp_stat(options, glob, params, 'Aprod', params.Aprod * 0.95, params.Aprod * 1.05, 5);
% cs_mu_bank            = comp_stat(options, glob, params, 'mu_bank', 0,0.5, 5);
% show_plots_cs(options, glob, params, cs_mu_bank)
% toc;

% load comp_stat.mat
% save comp_stat.mat


