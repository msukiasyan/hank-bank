function [paths, statst]  = do_transition_newton(opt, glob, p, init_state, final_ss)
    guesses0.Kt         = ones(p.Nt, 1) * final_ss.K;
    guesses0.x_at       = ones(p.Nt, 1) * final_ss.x_a;
    
    % Calculate the jacobian at the steady state (and at the provided paths
    % of exogenous variables)
    jac_for                 = zeros(2 * p.Nt, 2 * p.Nt);
    res_ss                  = transition_residuals(opt, glob, p, guesses0, init_state, final_ss);
    res_ss_vec              = [res_ss.K; res_ss.x_a];
    
    % wrt Kt
    for t = 1:p.Nt
        guesses0.Kt         = ones(p.Nt, 1) * final_ss.K;
        guesses0.Kt(t)      = guesses0.Kt(t) + opt.stepK_nt;
        
        res                 = transition_residuals(opt, glob, p, guesses0, init_state, final_ss);
        res_vec             = [res.K; res.x_a];
        jac_for(:, t)       = (res_vec - res_ss_vec) / opt.stepK_nt;
    end
    
    % wrt x_at
    for t = 1:p.Nt
        guesses0.x_at       = ones(p.Nt, 1) * final_ss.x_a;
        guesses0.x_at(t)    = guesses0.x_at(t) + opt.stepx_a_nt;
        
        res                 = transition_residuals(opt, glob, p, guesses0, init_state, final_ss);
        res_vec             = [res.K; res.x_a];
        jac_for(:, p.Nt + t)    = (res_vec - res_ss_vec) / opt.stepx_a_nt;
    end
    
    f                       = @(x) transition_residuals(opt, glob, p, x, init_state, final_ss);
    jac_for_inv             = inv(jac_for);

    trans_path              = [ones(p.Nt, 1) * final_ss.K; ones(p.Nt, 1) * final_ss.x_a];
    fval                    = f(struct('Kt', trans_path(1:p.Nt), 'x_at', trans_path(p.Nt+1:2*p.Nt)));
    fval_vec                = [fval.K; fval.x_a];
    for it = 1:opt.maxittrans_nt
        fnorm               = norm(fval_vec); 
 
        if opt.debug_trans >= 1
            disp(['iter ', int2str(it) ', norm = ' num2str(fnorm)  ]);
            figure(1); plot(p.tgrid, fval_vec(1:p.Nt)), ylim([-0.02 0.02]), grid; xlabel('Time (quarters)'); ylabel('Market clearing residual (K)') ;
        end

        if fnorm < opt.tol_nt
            break 
        end
        d                   = -(jac_for_inv * fval_vec);
        trans_path          = trans_path + d;
        fold_vec            = fval_vec;
        fval                = f(struct('Kt', trans_path(1:p.Nt), 'x_at', trans_path(p.Nt+1:2*p.Nt)));
        fval_vec            = [fval.K; fval.x_a];
        lu                  = jac_for_inv * (fval_vec - fold_vec);
        jac_for_inv         = jac_for_inv + ((d - lu) * (d' * jac_for_inv)) / (d' * lu);   
        
    end    
    paths.Kt                = trans_path(1:p.Nt);
    paths.x_at              = trans_path(p.Nt+1:2*p.Nt);
    [~, statst]             = f(struct('Kt', trans_path(1:p.Nt), 'x_at', trans_path(p.Nt+1:2*p.Nt)));
end