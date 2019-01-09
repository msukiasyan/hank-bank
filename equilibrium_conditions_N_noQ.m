% Inputs:  (1) vars: vector which contains the variables in the system, the time derivatives of 
%			         those variables, the expectational errors, and the shocks
%
% Outputs: (1) vResduals: residuals of equilibrium conditions, evaluated at vars

function vResidual = equilibrium_conditions_N_noQ(vars)

% global variables here
% Unpack vars
global alpha delta xi a1 theta f rho rhoB M ga
global Aswitch zzz aaa
global J zmin zmax amin amax I crit Delta the Var a z da dz aa zz mu s2 zmean sig2 Corr maxit P ownership mass
global varsSS rrhoTFP ssigmaTFP
global varsSS n_v n_g n_p n_shocks nEErrors nVars

V=vars(1:I*J*P) + varsSS(1:I*J*P);

g=vars(n_v+1:n_v+I*J*P-1) + varsSS(n_v+1:n_v+I*J*P-1);	% vector of distribution, removing last point		
g_end=1-sum(g);		% ensures that distribution integrates to 1
g = [g;g_end];
logAggregateTFP = vars(2*I*J*P);

% possibly need to rearrange stuff with dots
% what are the true state variables?
% e.g. is K a state variable or not?
% most likely alpha_hat and q_K are also jumpers with expectational errors

alpha_hat = vars(2*I*J*P+1) + varsSS(2*I*J*P+1);
N = vars(2*I*J*P+2) + varsSS(2*I*J*P+2);
K = vars(2*I*J*P+3) + varsSS(2*I*J*P+3);
r_plus = vars(2*I*J*P+4) + varsSS(2*I*J*P+4);


V = reshape(V,[I,J,P]);

VDot = vars(nVars+1:nVars+I*J*P);
gDot = vars(nVars+I*J*P+1:nVars+2*I*J*P-1);
logAggregateTFPDot = vars(nVars+2*I*J*P);

alpha_hat_dot =vars(nVars+2*I*J*P+1);
N_dot = vars(nVars+2*I*J*P+2);



VEErrors = vars(2*nVars+1:2*nVars+I*J*P);





alpha_hat_error = 0;

XN = f*N - M;
w = exp(logAggregateTFP) * (1 - alpha) * (K ^ alpha); 
r_minus = exp(logAggregateTFP) * alpha * K ^ (alpha -1) -delta; 

leverage = alpha_hat/theta;

lambda = leverage*(r_minus-r_plus);

% Initialize other variables, using vars to ensure everything is a dual number
residual_hjb = VEErrors;
gIntermediate = VEErrors;
gResidual = VEErrors(1:I*J*P-1);

Vaf = V;             
Vab = V;
Vzf = V;
Vzb = V;
Vzz = V;
c = V;



dividends = ownership*XN; %individual dividends
r = r_plus*(aaa>0) + r_minus*(aaa<0)+0;
%% ITERATION OF HJB

% forward difference
        Vaf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/da;
        Vaf(I,:,:) = (repmat(w*z,[1,1,P]) + r(I,:,:).*amax + dividends(I,:,:)).^(-ga); %will never be used, but impose state constraint a<=amax just in case
% backward difference
        Vab(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/da;
        Vab(1,:,:) = (repmat(w*z,[1,1,P]) + r(1,:,:).*amin +dividends(1,:,:)).^(-ga);  %state constraint boundary condition

        I_concave = Vab > Vaf;              %indicator whether value function is concave (problems arise if this is not the case)

        %consumption and savings with forward difference
        cf = Vaf.^(-1/ga);
        sf = w*zzz + r.*aaa - cf +dividends;
        %consumption and savings with backward difference
        cb = Vab.^(-1/ga);
        sb = w*zzz + r.*aaa - cb +dividends;
        %consumption and derivative of value function at steady state
        c0 = w*zzz + r.*aaa + dividends;
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
savings = c0 - c;

       %Now use a lopp
        %to iterate over different ownerships shares
        %almost surely it can be done without a loop
        %but I am lazy today
        A_big = [];
        for p = 1:P
        
        %CONSTRUCT MATRIX A
        X = - min(sb(:,:,p),0)/da;
        Y = - max(sf(:,:,p),0)/da + min(sb(:,:,p),0)/da;
        Z = max(sf(:,:,p),0)/da;
        
        updiag=[0]; %This is needed because of the peculiarity of spdiags.
        for j=1:J
            updiag=[updiag;Z(1:I-1,j);0];
        end
        
        centdiag=reshape(Y(:,:),I*J,1);
        
        lowdiag=X(2:I,1);
        for j=2:J
            lowdiag=[lowdiag;0;X(2:I,j)];
        end
        
        AA=spdiags(centdiag,0,I*J,I*J)+spdiags([updiag;0],1,I*J,I*J)+spdiags([lowdiag;0],-1,I*J,I*J);
        
        A = AA + Aswitch;
        %A_big = blkdiag(A_big,A);
        
        
        u_stacked = reshape(u(:,:,p),I*J,1);
        V_stacked = reshape(V(:,:,p),I*J,1);
        u_stacked + A * V_stacked;

        residual_hjb(((p-1)*I*J+1):p*I*J,1) =  u_stacked + A * V_stacked + VDot(((p-1)*I*J+1):p*I*J,1) + VEErrors(((p-1)*I*J+1):p*I*J,1) - rho * V_stacked;
        gIntermediate(((p-1)*I*J+1):p*I*J,1) = A' * g(((p-1)*I*J+1):p*I*J,1);
        if p == P
            gResidual(((p-1)*I*J+1):(p*I*J-1),1) = gDot(((p-1)*I*J+1):(p*I*J-1),1) - gIntermediate(((p-1)*I*J+1):(p*I*J-1),1);
        else
            gResidual(((p-1)*I*J+1):p*I*J,1) = gDot(((p-1)*I*J+1):p*I*J,1) - gIntermediate(((p-1)*I*J+1):p*I*J,1);
        end
        end

   

% HJB
%residual_hjb = reshape(u,I*J*P,1) + A_big * reshape(V,I*J*P,1) + VDot + VEErrors - rho * reshape(V,I*J*P,1);

% KFE 
% does it have to sum to 1? this has to be verified, because it seems it
% does not




reaa = repmat(reshape(aa,[I*J,1]),[P,1]);

B_plus = sum(g.*max(reaa,0));
B_minus = -sum(g.*min(reaa,0));

residual_alpha_hat = alpha_hat_error + (rhoB + f)*alpha_hat - f - alpha_hat*(r_plus+lambda) - alpha_hat_dot;
residual_N = N_dot - leverage*N*(r_minus-r_plus)-N*(r_plus-f) - M;
residual_assets = B_plus - (leverage-1)*N;
tfpResidual = logAggregateTFPDot + (1 - rrhoTFP) * logAggregateTFP - ssigmaTFP * 0;
vResidual=[residual_hjb;gResidual;tfpResidual;residual_alpha_hat;residual_N;residual_assets];








