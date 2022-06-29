
setPlotDefaults
%load sim_sfl/rs_H0.50_ni1e5wi1e5_mw1e+00_w1e-02.mat

fstr='sim_sfl/rs_H%.2f_ni1e5_wi%1.0e_mean.mat';
Hl = [0 1/4 1/2 .7 .71 .72 .73 .74 .75 1] 

figure
ha = tight_subplot(length(Hl)/2,2,[.02 .08],[.08 .03],[.11 .03]);
i2p = [1:2:length(Hl) 2:2:length(Hl)];
for i = 1:length(Hl)
    fns = sprintf(fstr,Hl(i),1e5); load(fns)
    axes(ha(i2p(i)))
    histogram(abs(cae),'BinWidth', 1/40, 'Normalization', 'pdf', 'EdgeColor','none','FaceColor',[.1 .1 .1]), hold on
    histogram(abs(ct),'BinWidth', 1/40, 'Normalization', 'pdf', 'EdgeColor','none', 'FaceColor',[.7 .7 .7])
    axis([0 2.6 0 2.5]) %axis([0 2.6 0 6000]) 
    if     Hl(i) == 0,    ylim([0 5])%ylim([0 1.1e4])
    elseif Hl(i) >= 0.75, ylim([0 11])  %ylim([0 3e4])
    end
    grid on, ca = gca; 
    if i<6, text(1.7, ca.YLim(2)*0.8, sprintf('H=%1.2f',Hl(i)), 'FontSize', 12)
    else,   text(0.1, ca.YLim(2)*0.8, sprintf('H=%1.2f',Hl(i)), 'FontSize', 12)
    end
    if i2p(i) < 9, xticklabels('')
    else,          xlabel('Activity')
    end
    if i == 3, ylabel('Probability density'), end
end


%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figAncsHists')
