%SET SOLUTION AND GRID PARAMETERS HERE
%This file creates all grids needed to solve the model

global J zmin zmax amin amax I crit Delta the Var a z da dz aa zz mu s2 zmean sig2 Corr maxit P aaa zzz
J=5;         % number of z points 
zmin = 0.7;   % Range z
zmax = 1.3;
amin = -0.1;    % borrowing constraint
amax = 30;    % range a
I=100;        % number of a points 
P = 1;       % number of points in the ownership distribution 

% ownership distribution
global mass ownership
mass = ones(1,P)/P; % specify mass of each point

% this section sets mass and ownership by hand
P = 1;
%mass=[0.99;0.01];
%ownership = [1;1];
P = 1;
mass = 1;
ownership = 1;
% this ends 
%ownership = linspace(0,5,P)./(sum(linspace(0,5,P).*mass)); %specify ownership distribution
ownership = reshape(ownership,[1,1,P]); %reshaping...
mass_vec = reshape(repmat(mass,I*J,1),[I*J*P,1]); 
mass = reshape(mass,[1,1,P]);
ownership = repmat(ownership,[I,J]);
mass = repmat(mass,[I,J]);


zmean = 1.0;      % mean O-U process (in levels). This parameter has to be adjusted to ensure that the mean of z (truncated gaussian) is 1.
sig2 = (0.10)^2;  % sigma^2 O-U
Corr = exp(-0.0);  % persistence -log(Corr)  O-U
%simulation parameters
maxit  = 100;     %maximum number of iterations in the HJB loop

crit = 10^(-10); %criterion HJB loop
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
aa = a*ones(1,J);
zz = ones(I,1)*z;
aaa = a.*ones(1,J,P);
zzz = ones(I,1,P).*z;
mu = the*(zmean - z);        %DRIFT (FROM ITO'S LEMMA)
s2 = sig2.*ones(1,J);        %VARIANCE (FROM ITO'S LEMMA)


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

global Aswitch
Aswitch=spdiags(centdiag,0,I*J,I*J)+spdiags(lowdiag,-I,I*J,I*J)+spdiags(updiag,I,I*J,I*J);


