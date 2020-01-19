function show_plots_cs(opt, glob, p, cs, stl)
    if nargin < 5
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
    title('Total Illiquid Assets');

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
    
    subplot(figure_height, figure_width, 12);
    plot(cs.cs_grid, cellfun(@(x) x.r_F, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Return');
    title('Return From Banks');
    
    subplot(figure_height, figure_width, 13);
    plot(cs.cs_grid, cellfun(@(x) x.w, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Wage');
    title('Wage');
    
    subplot(figure_height, figure_width, 14);
    plot(cs.cs_grid, cellfun(@(x) x.V_mean, cs.cs_stats), stl, 'LineWidth', 1.5)
    xlabel(cs.par);
    ylabel('Welfare');
    title('Welfare');
    
    if opt.GK && p.distGK == "twopoint"
        subplot(figure_height, figure_width, 15);
        yvar    = zeros(length(cs.cs_grid), 1);
        for cg  = 1:length(cs.cs_grid)
            yvar(cg)    = sum(cs.cs_stats{cg}.dtildea_vec .* ...
            cs.cs_stats{cg}.dtildeb_vec .* reshape(cs.cs_sol{cg}.dst, p.Nb * p.Na * p.Nz, 1) ...
            .* reshape(cs.cs_sol{cg}.cpol, p.Nb * p.Na * p.Nz, 1) .* (abs(cs.cs_stats{cg}.afrombaz - cs.cs_stats{cg}.a(2)) < 1e-9) ./ cs.cs_stats{cg}.fracGK);
        end
        plot(cs.cs_grid, yvar, stl, 'LineWidth', 1.5)
        xlabel(cs.par);
        ylabel('Welfare');
        title('Average Welfare of Bank Owners');
        
        subplot(figure_height, figure_width, 16);
        yvar    = zeros(length(cs.cs_grid), 1);
        for cg  = 1:length(cs.cs_grid)
            yvar(cg)    = sum(cs.cs_stats{cg}.dtildea_vec .* ...
            cs.cs_stats{cg}.dtildeb_vec .* reshape(cs.cs_sol{cg}.dst, p.Nb * p.Na * p.Nz, 1) ...
            .* reshape(cs.cs_sol{cg}.cpol, p.Nb * p.Na * p.Nz, 1) .* (abs(cs.cs_stats{cg}.afrombaz - cs.cs_stats{cg}.a(1)) < 1e-9) ./ (1 - cs.cs_stats{cg}.fracGK));
        end
        plot(cs.cs_grid, yvar, stl, 'LineWidth', 1.5)
        xlabel(cs.par);
        ylabel('Welfare');
        title('Average Welfare of the Rest');
    end
    
end