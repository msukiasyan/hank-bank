% This file finds stationary equilibrium of HANK with Families

clear all; close all; format long;
tic;

%--------------------------------------------------
%PARAMETERS
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

zmean = 1.0;      % mean O-U process (in levels). This parameter has to be adjusted to ensure that the mean of z (truncated gaussian) is 1.
sig2 = (0.10)^2;  % sigma^2 O-U
Corr = exp(-0.0);  % persistence -log(Corr)  O-U






relax = 0.99; % relaxation parameter 
J=20;         % number of z points 
zmin = 0.5;   % Range z
zmax = 1.5;
amin = 0;    % borrowing constraint
amax = 30;    % range a
I=100;        % number of a points 
P = 10;       % number of points in the ownership distribution

%simulation parameters
maxit  = 100;     %maximum number of iterations in the HJB loop
Delta = 1000;   %delta in HJB algorithm

%ORNSTEIN-UHLENBECK IN LEVELS
the = -log(Corr);
Var = sig2/(2*the);

%--------------------------------------------------
%VARIABLES 
a = linspace(amin,amax,I)';  %wealth vector
da = (amax-amin)/(I-1);      
z = linspace(zmin,zmax,J);   % productivity vector
dz = (zmax-zmin)/(J-1);
dz2 = dz^2;

mass = ones(1,P)/P;
ownership = linspace(0,3,P)./(sum(linspace(0,3,P).*mass));
ownership = reshape(ownership,[1,1,P]);
mass = reshape(mass,[1,1,P]);
ownership = repmat(ownership,[I,J]);
mass = repmat(mass,[I,J]);
aaa = a.*ones(1,J,P);
zzz = ones(I,1,P).*z;

mu = the*(zmean - z);        %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*ones(1,J,P);        %VARIANCE (FROM ITO'S LEMMA)

%Finite difference approximation of the partial derivatives
Vaf = zeros(I,J,P);             
Vab = zeros(I,J,P);
Vzf = zeros(I,J,P);
Vzb = zeros(I,J,P);
Vzz = zeros(I,J,P);
c = zeros(I,J,P);

%CONSTRUCT MATRIX Aswitch SUMMARIZING EVOLUTION OF z
yy = - s2/dz2 - mu/dz;
chi =  s2/(2*dz2);
zeta = mu/dz + s2/(2*dz2);

%This will be the upperdiagonal of the matrix Aswitch
updiag=zeros(I,1); %This is necessary because of the peculiar way spdiags is defined.
for j=1:J
    updiag=[updiag;repmat(zeta(j),I,1)];
end

%This will be the center diagonal of the matrix Aswitch
centdiag=repmat(chi(1)+yy(1),I,1);
for j=2:J-1
    centdiag=[centdiag;repmat(yy(j),I,1)];
end
centdiag=[centdiag;repmat(yy(J)+zeta(J),I,1)];

%This will be the lower diagonal of the matrix Aswitch
lowdiag=repmat(chi(2),I,1);
for j=3:J
    lowdiag=[lowdiag;repmat(chi(j),I,1)];
end

%Add up the upper, center, and lower diagonal into a sparse matrix
Aswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);
r_plus_max = f + M_N - 0.001;
iter = 1;
r_plus_vec = linspace(0.04529,0.045305,20);
%MAIN LOOP
for iter=1:20
    r_plus = r_plus_vec(iter);
disp('Main loop iteration')
disp(iter)

% SOLVE THE BANKING SECTOR

lambda = max(f - r_plus - M_N,0);
alpha_hat = f/(rhoB + f - r_plus - lambda);
r_minus = lambda*theta/alpha_hat + r_plus;
K = (alpha/(r_minus + delta))^(1/(1-alpha));
w = (1-alpha)*K^alpha;
leverage = alpha_hat/theta;
N = K/leverage;
XN = f*N - M_N*N;
v = log(w.*zz + XN + r_plus.*aa);


        % HAMILTON-JACOBI-BELLMAN EQUATION %
    for n=1:maxit
        V = v;
        % forward difference
        Vaf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
        Vaf(I,:) = (w*z + r_plus.*amax + XN).^(-ga); %will never be used, but impose state constraint a<=amax just in case
        % backward difference
        Vab(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
        Vab(1,:) = (w*z + r_plus.*amin + XN).^(-ga);  %state constraint boundary condition

        I_concave = Vab > Vaf;              %indicator whether value function is concave (problems arise if this is not the case)

        %consumption and savings with forward difference
        cf = Vaf.^(-1/ga);
        sf = w*zz + r_plus.*aa - cf + XN;
        %consumption and savings with backward difference
        cb = Vab.^(-1/ga);
        sb = w*zz + r_plus.*aa - cb + XN;
        %consumption and derivative of value function at steady state
        c0 = w*zz + r_plus.*aa + XN;
        Va0 = c0.^(-ga);

        % dV_upwind makes a choice of forward or backward differences based on
        % the sign of the drift
        If = sf > 0; %positive drift --> forward difference
        Ib = sb < 0; %negative drift --> backward difference
        I0 = (1-If-Ib); %at steady state
        %make sure backward difference is used at amax
        %     Ib(I,:) = 1; If(I,:) = 0;
        %STATE CONSTRAINT at amin: USE BOUNDARY CONDITION UNLESS sf > 0:
        %already taken care of automatically

        Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term

        c = Va_Upwind.^(-1);
        u = log(c);

        %CONSTRUCT MATRIX A
        X = - min(sb,0)/da;
        Y = - max(sf,0)/da + min(sb,0)/da;
        Z = max(sf,0)/da;
        
        updiag=[0]; %This is needed because of the peculiarity of spdiags.
        for j=1:J
            updiag=[updiag;Z(1:I-1,j);0];
        end
        
        centdiag=reshape(Y,I*J,1);
        
        lowdiag=X(2:I,1);
        for j=2:J
            lowdiag=[lowdiag;0;X(2:I,j)];
        end
        
        AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
        
        A = AA + Aswitch;
        B = (1/Delta + rho)*speye(I*J) - A;

        u_stacked = reshape(u,I*J,1);
        V_stacked = reshape(V,I*J,1);

        b = u_stacked + V_stacked/Delta;

        V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS

        V = reshape(V_stacked,I,J);

        Vchange = V - v;
        v = V;

        dist(n) = max(max(abs(Vchange)));
        if dist(n)<crit
            disp('Value Function Converged, Iteration = ')
            disp(n)
            break
        end
    end
    toc;
    % FOKKER-PLANCK EQUATION %
    AT = A';
    b = zeros(I*J,1);

    %need to fix one value, otherwise matrix is singular
    i_fix = 1;
    b(i_fix)=.1;
    row = [zeros(1,i_fix-1),1,zeros(1,I*J-i_fix)];
    AT(i_fix,:) = row;

    %Solve linear system
    gg = AT\b;
    g_sum = gg'*ones(I*J,1)*da*dz;
    gg = gg./g_sum;

    g = reshape(gg,I,J);
    
    % Update aggregate capital
    B_plus = sum(g'*a*da*dz);
    disp(B_plus)
   
err(iter) = (B_plus - (leverage-1)*N);

end
%GRAPHS 
% 
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
