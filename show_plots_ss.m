function show_plots_ss(opt, glob, p, st, stl)
    if nargin < 4
        stl = '-';
    end
  
    %% Plot distributions
    
    figure(1);
    sgtitle('Distributions');
    
    for nz = 1:p.Nz
        subplot(p.Nz, 2, nz * 2 - 1);
        plot(p.b, st.ldist{nz}, stl, 'LineWidth', 1.5);
        xlabel('Liquid Wealth');
        ylabel('Density');
        title(['Liquid Distribution, Type ' num2str(nz)]);

        subplot(p.Nz, 2, nz * 2);
        plot(p.a(1:27), st.ildist{nz}(1:27), stl, 'LineWidth', 1.5);
        xlabel('Illiquid Wealth');
        ylabel('Density');
        title(['Illiquid Distribution, Type ' num2str(nz)]);
    end

%     subplot(2,2,3);
%     plot(p.b, st.ldist{2}, stl, 'LineWidth', 1.5);
%     xlabel('Liquid Wealth');
%     ylabel('Density');
%     title('Liquid Distribution, High Type');
% 
%     subplot(2,2,4);
%     plot(p.a, st.ildist{2}, stl, 'LineWidth', 1.5);
%     xlabel('Illiquid Wealth');
%     ylabel('Density');
%     title('Illiquid Distribution, High Type');
end