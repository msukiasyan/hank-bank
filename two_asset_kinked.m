clear all; close all; clc;
global chi0 chi1 ga rho xi rb_spread Aprod alp delt I J crit da db a b Nz bmin amax z la_mat w
global Delta maxit bmax amin
global rho_bank f_bank theta_bank
ga = 2; %CRRA utility with parameter gamma
rho = 0.06; %discount rate
chi0 = 0.03;%0.01;
chi1 = 2.0;
xi = 0.1; %fraction of income that is automatically deposited
rb_spread = 0.025;
Aprod = 1.0;
alp = 0.4;
delt = 0.05;

rho_bank    = -log(0.96);
f_bank      = 0.05;
theta_bank  = 0.19;

%Income process (two-state Poisson process):
w = 4;
Nz = 2;
z      = [.8,1.3];
%z      = [.7,1.4];
la_mat = [-1/3, 1/3; 1/3, -1/3];

crit = 10^(-5);
Delta = 100;
maxit = 35;

%grids
I = 500;
bmin = -2;
%bmin = 0;
bmax = 47.9;
b = linspace(bmin,bmax,I)';
db = (bmax-bmin)/(I-1);

J= 50;
amin = 0;
amax = 100;
a = linspace(amin,amax,J);
da = (amax-amin)/(J-1);

[cpol, dpol, bpol, apol, dst] = get_policies(0.08, 0.06, 0.0575);
TL = 0;
TK = 0;
ldist = zeros(Nz, 1);
ldist = sum(dst(:, :, 1), 2) * da;
plot(b, ldist)
return

for nz = 1:Nz
    TL = TL + sum(dst(:, :, nz), 2)' * b * db * da;
end

for nz = 1:Nz
    TK = TK + sum(dst(:, :, nz), 1) * a' * da * db;
end

TL
TK

sol = fsolve(@eqs, [11, 0.045], optimoptions('fsolve', 'Display', 'iter'))

%%%%%%%%%%%
% FIGURES %
%%%%%%%%%%%


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


