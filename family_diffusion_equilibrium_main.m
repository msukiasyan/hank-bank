% This file finds stationary equilibrium of HANK with Families

clear all; close all; format long;
tic;

%--------------------------------------------------
%
set_parameters;
set_grids;

global r_plus_ss kappa_ss r_minus_ss K_ss N_ss B_plus_ss B_minus_ss alpha_hat_ss...
        lambda_ss XN_ss V_ss g_ss leverage_ss 
out = fsolve(@steady_state,[0.03;-2])

% %
% % %SAVINGS POLICY FUNCTION
% figure
% ss = w*zz + r.*aa - c;
% icut = 50;
% acut = a(1:icut);
% sscut = ss(1:icut,:);
% set(gca,'FontSize',14)
% surf(acut,z,sscut')
% view([45 25])
% xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
% ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
% zlabel('Savings $s(a,z)$','FontSize',14,'interpreter','latex')
% xlim([amin max(acut)])
% ylim([zmin zmax])
% 
% %WEALTH DISTRIBUTION
% figure
% icut = 50;
% acut = a(1:icut);
% gcut = g(1:icut,:);
% set(gca,'FontSize',14)
% surf(acut,z,gcut')
% view([45 25])
% xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
% ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
% zlabel('Density $f(a,z)$','FontSize',14,'interpreter','latex')
% xlim([amin max(acut)])
% ylim([zmin zmax])
