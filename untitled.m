% return residuals of equilibrium conditions
function [outputArg1,outputArg2] = steady_state(inp)

global alpha delta xi a1 theta f rho rhoB M_N ga
global Aswitch
global J zmin zmax amin amax I crit Delta the Var a z da dz aa zz mu s2
% load inputs
r_plus = inp(1);
kappa = inp(2);



%% Aggregate equations %%

lambda = max(f - r_plus - M_N,0); %multiplier on the budget constraint 
alpha_hat = f/(rhoB + f - r_plus - lambda); %marginal value of net worth
r_minus = lambda*theta/alpha_hat + r_plus; %interest on loans
K = (alpha/(r_minus + delta))^(1/(1-alpha)); %capital
w = (1-alpha)*K^alpha; %real wage
leverage = alpha_hat/theta; %leverage
N = K/(kappa*leverage); %recover net worth 
XN = f*N - M_N*N; %aggregate dividend
dividends = ownership*XN; %individual dividends
v = log(w.*zzz + dividends + r_plus.*aaa);
r = r_plus*(aaa>=0) + r_minus*(aaa<0);

%%  HAMILTON-JACOBI-BELLMAN EQUATION %%
    for n=1:maxit
        V = v;
        % forward difference
        Vaf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/da;
        Vaf(I,:,:) = (w*z + r(I,:,:).*amax + dividends(I,:,:)).^(-ga); %will never be used, but impose state constraint a<=amax just in case
        % backward difference
        Vab(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))/da;
        Vab(1,:,:) = (w*z + r(I,:,:).*amin +dividends(I,:,:)).^(-ga);  %state constraint boundary condition

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
        %Now use a lopp
        %to iterate over different ownerships shares
        %almost surely it can be done without a loop
        %but I am lazy today
        AT = [];
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
        B = (1/Delta + rho)*speye(I*J) - A;

        u_stacked = reshape(u(:,:,p),I*J,1);
        V_stacked = reshape(V(:,:,p),I*J,1);

        b = u_stacked + V_stacked/Delta;

        V_stacked = B\b; %SOLVE SYSTEM OF EQUATIONS

        V(:,:,p) = reshape(V_stacked,I,J);
                AT = blkdiag(AT,A');
        end
        
        Vchange = V - v;
        v = V;

        dist(n) = max(max(max(abs(Vchange))));

        if dist(n)<crit
            disp('Value Function Converged, Iteration = ')
            disp(n)
            break
        end

    toc;
    end
    for p = 1:P
   % b = zeros(I*J*P);
    ATemp = AT(I*J*(p-1)+1:I*J*p,I*J*(p-1)+1:I*J*p);
    %need to fix one value, otherwise matrix is singular
     b = zeros(I*J,1);

    %need to fix one value, otherwise matrix is singular
    i_fix = 1;
    b(i_fix)=.1;
    row = [zeros(1,i_fix-1),1,zeros(1,I*J-i_fix)];
    ATemp(i_fix,:) = row;


    %Solve linear system
    gg = ATemp\b;
    g_sum = gg'*ones(I*J,1)*da*dz*mass(p);
    gg = gg./g_sum;
    g(:,:,p) = reshape(gg,I,J);
    end

    
    % Update aggregate capital
    B_plus =  sum(sum(sum(g.*aaa*da*dz.*mass)));
    disp(B_plus)