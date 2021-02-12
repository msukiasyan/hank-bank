function show_plots_distributional_irf(opt, glob, p, stats, paths, statst)
    
    figure;
    % sgtitle('1% capital stock shock, AR factor 0.95');
    %sgtitle('5% net worth shock');
    sgtitle('Negative Taylor rule innovation');

    subplot(2, 5, 1);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.cons_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Consumption');
    
    
    subplot(2, 5, 2);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.dep_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference');
    title('Deposits');
    
    
    subplot(2, 5, 3);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.h_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Labor');

    subplot(2, 5, 4);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.disposable_inc_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Disposable income');

        subplot(2, 5, 5);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.s_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('S');
    
        subplot(2, 5, 6);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.m_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('M');
    
         subplot(2, 5, 7);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.liq_inc_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('Liquid asset income');
    
         subplot(2, 5, 8);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.illiq_inc_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('Illiquid asset income');
    
         subplot(2, 5, 9);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.b_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('Liquid assets');
    
             subplot(2, 5, 10);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.illiq.a_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('Illiquid assets');
    
    
    set(findobj(gcf,'type','axes'), 'XGrid', 'on', 'YGrid', 'on');
    
    fig = gcf;
    fig.Position(3) = fig.Position(3) + 250;
% add legend
    Lgnd = legend({'0th','10th', '25th', '50th', '75th', '90th','95th'});
    Lgnd.Position(1) = 0.01;
    Lgnd.Position(2) = 0.4;

%% now liquid
    figure;
    % sgtitle('1% capital stock shock, AR factor 0.95');
    %sgtitle('5% net worth shock');
    sgtitle('Negative Taylor rule innovation');

    subplot(2, 5, 1);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.cons_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Consumption');
    
    
    subplot(2, 5, 2);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.dep_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference');
    title('Deposits');
    
    
    subplot(2, 5, 3);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.h_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Labor');

    subplot(2, 5, 4);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.disposable_inc_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Percent dev.');
    title('Disposable income');

        subplot(2, 5, 5);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.s_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('S');
    
        subplot(2, 5, 6);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.m_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('M');
    
         subplot(2, 5, 7);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.liq_inc_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('Liquid asset income');
    
         subplot(2, 5, 8);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.illiq_inc_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('Illiquid asset income');
    
         subplot(2, 5, 9);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.b_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('Liquid assets');
    
             subplot(2, 5, 10);
    for k  = 1:7
    plot(p.tgrid(1:26), (cellfun(@(x) x.liq.a_irf_percentile(k), statst(1:26)) ), 'LineWidth', 1.5);
    hold on
    end
    xlabel('Quarters');
    ylabel('Difference.');
    title('Illiquid assets');
    
    
    set(findobj(gcf,'type','axes'), 'XGrid', 'on', 'YGrid', 'on');
    
    fig = gcf;
    fig.Position(3) = fig.Position(3) + 250;
% add legend
    Lgnd = legend({'0th','10th', '25th', '50th', '75th', '90th','95th'});
    Lgnd.Position(1) = 0.01;
    Lgnd.Position(2) = 0.4;

end