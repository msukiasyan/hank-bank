function show_plots_cs(opt, p, cs, stl)
    if nargin < 4
        stl = '-';
    end
  
    %% Plot comparative statics
    figure_width    = 4;
    figure_height   = 4;
    
    figure;
    sgtitle(['Comparative Statics wrt ' cs.par]);
    
    subplot(figure_height, figure_width, 1);
    plot(cs.cs_grid, cellfun(@(x) x.cons_mean, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Total Consumption');
    title('Total Consumption');

    subplot(figure_height, figure_width, 2);
    plot(cs.cs_grid, cellfun(@(x) x.cons_gini, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Gini of Consumption');
    title('Gini of Consumption');

    subplot(figure_height, figure_width, 3);
    plot(cs.cs_grid, cellfun(@(x) x.a_mean, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Total Illiquid Assets');
    title('Total Illiquid Assets = NW');

    subplot(figure_height, figure_width, 4);
    plot(cs.cs_grid, cellfun(@(x) x.a_gini, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Gini of Illiquid Assets');
    title('Gini of Illiquid Assets');
    
    subplot(figure_height, figure_width, 5);
    plot(cs.cs_grid, cellfun(@(x) x.b_mean, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Total Liquid Assets');
    title('Total Liquid Assets');
    
    subplot(figure_height, figure_width, 6);
    plot(cs.cs_grid, cellfun(@(x) x.b_gini, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Gini of Liquid Assets');
    title('Gini of Liquid Assets');
    
    subplot(figure_height, figure_width, 7);
    plot(cs.cs_grid, cellfun(@(x) x.K, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Capital');
    title('Capital');
    
    subplot(figure_height, figure_width, 8);
    plot(cs.cs_grid, cellfun(@(x) x.x_a, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Leverage');
    title('Leverage');
    
    subplot(figure_height, figure_width, 9);
    plot(cs.cs_grid, cellfun(@(x) x.spread, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Spread');
    title('Spread');
    
    subplot(figure_height, figure_width, 10);
    plot(cs.cs_grid, cellfun(@(x) x.r_plus, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Deposit Rate');
    title('Deposit Rate');
    
    subplot(figure_height, figure_width, 11);
    plot(cs.cs_grid, cellfun(@(x) x.r_minus, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Lending Rate');
    title('Lending Rate');
    
end