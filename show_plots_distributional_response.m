function show_plots_distributional_response(opt, glob, p, stats, paths, statst,distr_response)
    
    figure;
  
    sgtitle('Negative Taylor rule innovation');

    subplot(2, 5, 1);
    plot(linspace(0,0.95,20), distr_response.illiq.cons, 'LineWidth', 1.5);
    hold on

    xlabel('Percentile');
    ylabel('Percent dev.');
    title('Consumption');
    
    
    subplot(2, 5, 2);
    plot(linspace(0,0.95,20), distr_response.illiq.h, 'LineWidth', 1.5);
    hold on
 
    xlabel('Percentile');
    ylabel('Percent dev.');
    title('Labor');
    
    
    subplot(2, 5, 3);
    plot(linspace(0,0.95,20), distr_response.illiq.d, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('Illiqud Deposit');
    
    subplot(2, 5, 4);
    plot(linspace(0,0.95,20), distr_response.illiq.m, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('M');
    
    
    subplot(2, 5, 5);
    plot(linspace(0,0.95,20), distr_response.illiq.s, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('S');
    
        subplot(2, 5, 6);
    plot(linspace(0,0.95,20), distr_response.illiq.labor_inc, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Percent dev.');
    title('Labor income');
    
        subplot(2, 5, 7);
    plot(linspace(0,0.95,20), distr_response.illiq.liq_inc, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('Liquid income');
    
    subplot(2, 5, 8);
    plot(linspace(0,0.95,20), distr_response.illiq.illiq_inc, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('Illiquid income');

        subplot(2, 5, 9);
    plot(linspace(0,0.95,20), distr_response.illiq.disposable_inc, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Percent dev.');
    title('Disposable income');
    
       subplot(2, 5, 10);
    plot(linspace(0,0.95,20), distr_response.illiq.liquid, 'LineWidth', 1.5,'Color','r');
    hold on
    
    xlabel('Percentile');
    ylabel('Level');
    title('Liquid assets');
    

    
    %% LIQUID
     figure;
  
    sgtitle('Negative Taylor rule innovation');

    subplot(2, 5, 1);
    plot(linspace(0,0.95,20), distr_response.liq.cons, 'LineWidth', 1.5);
    hold on

    xlabel('Percentile');
    ylabel('Percent dev.');
    title('Consumption');
    
    
    subplot(2, 5, 2);
    plot(linspace(0,0.95,20), distr_response.liq.h, 'LineWidth', 1.5);
    hold on
 
    xlabel('Percentile');
    ylabel('Percent dev.');
    title('Labor');
    
    
    subplot(2, 5, 3);
    plot(linspace(0,0.95,20), distr_response.liq.d, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('Illiqud Deposit');
    
    subplot(2, 5, 4);
    plot(linspace(0,0.95,20), distr_response.liq.m, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('M');
    
    
    subplot(2, 5, 5);
    plot(linspace(0,0.95,20), distr_response.liq.s, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('S');
    
        subplot(2, 5, 6);
    plot(linspace(0,0.95,20), distr_response.liq.labor_inc, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Percent dev.');
    title('Labor income');
    
        subplot(2, 5, 7);
    plot(linspace(0,0.95,20), distr_response.liq.liq_inc, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('Liquid income');
    
    subplot(2, 5, 8);
    plot(linspace(0,0.95,20), distr_response.liq.illiq_inc, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Difference');
    title('Illiquid income');

        subplot(2, 5, 9);
    plot(linspace(0,0.95,20), distr_response.liq.disposable_inc, 'LineWidth', 1.5);
    hold on
    
    xlabel('Percentile');
    ylabel('Percent dev.');
    title('Disposable income');
    
    
       subplot(2, 5, 10);
    plot(linspace(0,0.95,20), distr_response.liq.illiquid, 'LineWidth', 1.5,'Color','r');
    hold on
    
    xlabel('Percentile');
    ylabel('Level');
    title('Illiquid assets');
end