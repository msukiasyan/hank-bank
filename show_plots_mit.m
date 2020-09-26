function show_plots_mit(opt, glob, p, stats, paths, statst)
    
    figure;
    % sgtitle('1% capital stock shock, AR factor 0.95');
    sgtitle('1% productivity shock, \rho = 0.95');

    subplot(4, 4, 1);
    plot(p.tgrid(1:26), cellfun(@(x) x.I, statst(1:26)) / (stats.K * p.delta) * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Investment');
    
    subplot(4, 4, 2);
    plot(p.tgrid(1:26), cellfun(@(x) x.Y, statst(1:26)) / statst{p.Nt}.Y * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Output');
    
    subplot(4, 4, 3);
    plot(p.tgrid(1:26), cellfun(@(x) x.cons_mean, statst(1:26)) / stats.cons_mean * 100 - 100, ...
        'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Consumption');
    
    subplot(4, 4, 4);
    plot(p.tgrid(1:26), cellfun(@(x) x.N, statst(1:26)) / stats.N * 100 - 100, ...
        'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Agg Labor');
    
    subplot(4, 4, 5);
    plot(p.tgrid(1:26), cellfun(@(x) x.TD, statst(1:26)) / stats.TD * 100 - 100, ...
        'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Liquid Assets');
    
    subplot(4, 4, 6);
    plot(p.tgrid(1:26), cellfun(@(x) x.TS, statst(1:26)) / stats.TS * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Illiquid assets');
    
    subplot(4, 4, 7);
    plot(p.tgrid(1:26), (cellfun(@(x) x.r_F, statst(1:26)) - stats.r_F) * 400, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev. (annualized)');
    title('Return on bank equity');
    
    subplot(4, 4, 8);
    plot(p.tgrid(1:26), cellfun(@(x) x.NW, statst(1:26)) / stats.NW * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Net Worth');
    
    subplot(4, 4, 9);
    plot(p.tgrid(1:26), cellfun(@(x) x.TD_bank, statst(1:26)) / stats.TD_bank * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Total deposits');

    subplot(4, 4, 10);
    plot(p.tgrid(1:26), paths.x_at(1:26), 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Level');
    title('Leverage');

    subplot(4, 4, 11);
    plot(p.tgrid(1:26), (cellfun(@(x) x.r_minus, statst(1:26)) - stats.r_minus) * 400, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev. (annualized)');
    title('Borrowing rate');

    subplot(4, 4, 12);
    plot(p.tgrid(1:26), (cellfun(@(x) x.r_plus, statst(1:26)) - stats.r_plus) * 400, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev. (annualized)');
    title('Deposit (liquid) rate');

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
    
    subplot(4, 4, 13);
    plot(p.tgrid(1:26), (cellfun(@(x) x.spread, statst(1:26)) - stats.spread) * 400, ...
        'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev. (annualized)');
    title('Spread');
    
    subplot(4, 4, 14);
    plot(p.tgrid(1:26), cellfun(@(x) x.q, statst(1:26)) / 1 * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Capital price');
    
    subplot(4, 4, 15);
    plot(p.tgrid(1:26), cellfun(@(x) x.cons_gini, statst(1:26)) / stats.cons_gini * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Consumption Gini');
    
    subplot(4, 4, 16);
    plot(p.tgrid(1:26), cellfun(@(x) x.c90 - x.c10, statst(1:26)) / (stats.c90 - stats.c10) * 100 - 100, 'LineWidth', 1.5);
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Consumption 90-10 percentile');
    
    set(findobj(gcf,'type','axes'), 'XGrid', 'on', 'YGrid', 'on');

end