function show_plots_mit(opt, glob, p, stats, paths, statst)
    
    figure;
    % sgtitle('1% capital stock shock, AR factor 0.95');
    sgtitle('5% net worth shock');

    % subplot(4, 3, 1);
    % plot(p.tgrid(1:26), -0.01 * exp(-0.05 * p.tgrid(1:26)) * 100, 'LineWidth', 1.5);
    % xlabel('Quarters');
    % ylabel('Percent dev.');
    % title('Productivity');
%     subplot(4, 3, 1);
%     plot(p.tgrid(1:26), cellfun(@(x) x.V_mean, statst(1:26)) / stats.V_mean * 100 - 100, ...
%         p.tgrid(1:26), cellfun(@(x) x.V_mean_p, statst(1:26)) / stats.V_mean_p * 100 - 100, ...
%         p.tgrid(1:26), cellfun(@(x) x.V_mean_w, statst(1:26)) / stats.V_mean_w * 100 - 100, 'LineWidth', 1.5);
%     legend('Mean', 'Poor', 'Wealthy');
%     xlabel('Quarters');
%     ylabel('Percent dev.');
%     title('Value');

    subplot(4, 4, 1);
    plot(p.tgrid(1:26), paths.Kt(1:26) / stats.K * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Capital');
    
    subplot(4,4, 2);
    plot(p.tgrid(1:26), cellfun(@(x) x.Y, statst(1:26)) / statst{p.Nt}.Y * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Output');

    subplot(4, 4, 3);
    plot(p.tgrid(1:26), paths.x_at(1:26) / stats.x_a * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Leverage');

    subplot(4, 4, 4);
    plot(p.tgrid(1:26), (cellfun(@(x) x.r_minus, statst(1:26)) - stats.r_minus) * 400, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev. (annualized)');
    title('Borrowing rate');

    subplot(4, 4, 5);
    plot(p.tgrid(1:26), (cellfun(@(x) x.r_plus, statst(1:26)) - stats.r_plus) * 400, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev. (annualized)');
    % ylim([-10 50]);
    title('Deposit (liquid) rate');

%     subplot(4, 3, 5);
%     plot(p.tgrid(1:26), ((paths.x_at(1:26) .* cellfun(@(x) x.spread, statst(1:26)) ...
%                                 + cellfun(@(x) x.r_plus, statst(1:26))) - stats.r_F) * 400, 'LineWidth', 1.5);
%     xlabel('Quarters');
%     ylabel('Percent dev. (annualized)');
%     % ylim([-30 40]);
%     title('Bank (illiquid) return');

    subplot(4, 4, 6);
    plot(p.tgrid(1:26), cellfun(@(x) x.cons_mean, statst(1:26)) / stats.cons_mean * 100 - 100, ...
        'LineWidth', 1.5);
        %p.tgrid(1:26), cellfun(@(x) x.cons_mean_p, statst(1:26)) / stats.cons_mean_p * 100 - 100, ...
        %p.tgrid(1:26), cellfun(@(x) x.cons_mean_w, statst(1:26)) / stats.cons_mean_w * 100 - 100, ...
        
    % legend('Mean', 'Poor', 'Wealthy');
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Consumption');

    subplot(4, 4, 7);
    plot(p.tgrid(1:26), cellfun(@(x) x.TD, statst(1:26)) / stats.TD * 100 - 100, ...
        'LineWidth', 1.5);
        %p.tgrid(1:26), cellfun(@(x) x.b_mean_p, statst(1:26)) / stats.b_mean_p * 100 - 100, ...
        %p.tgrid(1:26), cellfun(@(x) x.b_mean_w, statst(1:26)) / stats.b_mean_w * 100 - 100, ...
        
    % plot(p.tgrid(1:26), cellfun(@(x) x.TD, statst(1:26)) / stats.TD * 100 - 100, 'LineWidth', 1.5);
    % legend('Mean', 'Poor', 'Wealthy');
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Liquid Assets');

    subplot(4, 4, 8);
    plot(p.tgrid(1:26), cellfun(@(x) x.TS, statst(1:26)) / stats.TS * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Net Worth');
    
%     figure;
%     sgtitle('1% capital stock shock, inequality');

    subplot(4, 4, 9);
    plot(p.tgrid(1:26), cellfun(@(x) x.cons_gini, statst(1:26)) / stats.cons_gini * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Consumption Gini');
    
    subplot(4, 4, 10);
    plot(p.tgrid(1:26), cellfun(@(x) sqrt(x.cons_var), statst(1:26)) / sqrt(stats.cons_var) * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Consumption dispersion (std dev)');

%     subplot(4, 3, 3);
%     plot(p.tgrid(1:26), cellfun(@(x) x.b_gini, statst(1:26)) / stats.b_gini * 100 - 100, 'LineWidth', 1.5);
%     xlabel('Quarters');
%     ylabel('Percent dev.');
%     title('Liquid Gini');
% 
%     subplot(4, 3, 4);
%     plot(p.tgrid(1:26), cellfun(@(x) x.a_gini, statst(1:26)) / stats.a_gini * 100 - 100, 'LineWidth', 1.5);
%     xlabel('Quarters');
%     ylabel('Percent dev.');
%     title('Illiquid Gini');
    
    subplot(4, 4, 11);
    plot(p.tgrid(1:26), cellfun(@(x) x.total_inc_mean, statst(1:26)) / stats.total_inc_mean * 100 - 100, ...
        'LineWidth', 1.5);
        % p.tgrid(1:26), cellfun(@(x) x.total_inc_mean_p, statst(1:26)) / stats.total_inc_mean_p * 100 - 100, ...
        % p.tgrid(1:26), cellfun(@(x) x.total_inc_mean_w, statst(1:26)) / stats.total_inc_mean_w * 100 - 100, ...
        
    % legend('Mean', 'Poor', 'Wealthy');
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Total income');
    
    subplot(4, 4, 12);
    plot(p.tgrid(1:26), cellfun(@(x) x.I, statst(1:26)) / (p.delta * stats.K) * 100 - 100, ...
        'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Investment');
    
        subplot(4, 4, 13);
    plot(p.tgrid(1:26), cellfun(@(x) x.w, statst(1:26)) / (stats.w) * 100 - 100, ...
        'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Wage');
    
        
        subplot(4, 4, 14);
    plot(p.tgrid(1:26), cellfun(@(x) x.piP, statst(1:26))* 400 , ...
        'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Annualized Rate.');
    title('Inflation');
    
            subplot(4, 4, 14);
    plot(p.tgrid(1:26), (cellfun(@(x) x.r_plusnom, statst(1:26)) - stats.r_plus) * 400, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Annualized Rate.');
    title('Nominal interest rate');
    
    
    set(findobj(gcf,'type','axes'), 'XGrid', 'on', 'YGrid', 'on');

end