% This file finds stationary equilibrium of HANK with Families

clear all; close all; format long;
tic;

%--------------------------------------------------
%
set_parameters;
set_grids;
global r_plus_ss kappa_ss r_minus_ss K_ss N_ss B_plus_ss B_minus_ss alpha_hat_ss...
        lambda_ss XN_ss V_ss g_ss leverage_ss 
global varsSS n_v n_g n_p n_shocks nEErrors nVars

n_v = I*J*P;              % number of jump variables
n_g = I*J*P;              % number of endogeneous state variables
n_p = 4;              % number of static relations 
n_shocks = 1;           % 
nEErrors = n_v+2;
nVars = n_v + n_g + n_p;




out = fsolve(@steady_state_N,[0.03;-2;1])

% prepare steady state variables
varsSS(1:I*J*P) = reshape(V_ss,[I*J*P,1]);
varsSS(n_v+1:n_v+I*J*P) = reshape(g_ss,[I*J*P,1]);
varsSS(2*I*J*P+1) = alpha_hat_ss;
varsSS(2*I*J*P+2) = N_ss;
varsSS(2*I*J*P+3) = K_ss;
varsSS(2*I*J*P+4) = r_plus_ss;
varsSS = reshape(varsSS,[nVars,1]);


%% Step 2: Linearize Model Equations

fprintf('Taking derivatives of equilibrium conditions...\n');
t0 = tic;


% Prepare automatic differentiation
vars = zeros(nVars + nVars + nEErrors + n_shocks,1);
vars = myAD(vars);
resids = equilibrium_conditions_N_noQ(vars);



% Evaluate derivatives

% Extract derivative values
derivs = getderivs(resids);

t_derivs = toc(t0);
fprintf('Time to compute derivatives: %2.4f seconds\n\n\n',t_derivs);
if t_derivs > 1
    warning('If you do not compile mex/C files for the automatic differentiation, matrix vector multiplication will be slow');
    disp('Press any key to continue...');
    pause();
end

% %% Step 3: Solve out Static Constraints or Reduce the Model
% % Extract derivatives
g1 = -derivs(:,1:nVars);
g0 = derivs(:,nVars+1:2*nVars);
pi = -derivs(:,2*nVars+1:2*nVars+nEErrors);
psi = -derivs(:,2*nVars+nEErrors+1:2*nVars+nEErrors+n_shocks);
constant = sparse(nVars,1);
% 
% % State Variables
%     % Solve out static constraints
    fprintf('Solving Out Static Constraints ...\n');
    [state_red,inv_state_red,g0,g1,constant,pi,psi] = clean_G0_sparse(g0,g1,constant,pi,1);
    n_g_red = n_g;
% 
%     % Create identity matrix for code reuse below
%     from_spline = speye(n_g_red + n_v);
%     to_spline = speye(n_g_red + n_v);
%     n_splined = n_v;
% 
% 
% 
% %% Step 4: Solve Linear Systems
% t0 = tic;
% fprintf('Solving linear system...\n');
% [G1, ~, impact, eu, F] = schur_solver(g0,g1,constant,psi,pi,1,1,1,2);
% fprintf('...Done!\n')
% fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
% fprintf('Time to solve linear system: %2.4f seconds\n\n\n',toc(t0));