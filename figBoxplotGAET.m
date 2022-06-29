
lw=1.5;

figure, ha = tight_subplot(2,1,[.08],[.08 .03],[.11 .03]);
axes(ha(1)) %subplot(211)
xlabel('Level'), ylabel('Activity / Gain')
axis([.8 7.2 0 1]), grid on, hold on, box on
ph(1)=errorbar(mnaA./sdE,semaA./sdE,'--', 'LineWidth', lw);
ph(2)=errorbar(sdA./sdE,sesdA./sdE,'--', 'LineWidth', lw);
ph(3)=plot(Pfgi,'--', 'LineWidth', lw);
ph(4)=errorbar(mnT./sdE,semaT./sdE,'--', 'LineWidth', lw);
ph(5)=errorbar(mnaE./sdE,semaE./sdE,'--', 'LineWidth', lw);
ph(6:10)=line([1;7]*ones(1,5), [1;1]*[maa sda gst tht mae], 'LineWidth', 1);
legend(ph(1:10), {'$\langle|A_l|\rangle$','$\sigma_{A_l}$','g','$\langle\Theta_l\rangle$','$\langle|\mathcal{E}_l|\rangle$',...
    '$\langle|A_*|\rangle$','$\sigma_{A_*}$','$g_*$','$\Theta_*$','$\langle|\mathcal{E}_*|\rangle$'}, ...
    'Interpreter', 'latex', 'Location','eastoutside')
axes(ha(2)) %subplot(212)
xlabel('Level'), ylabel('Normalized activity / gain')
xlim([.8 7.2]), grid on, hold on, box on
errorbar(mnaA./sdE/maa,semaA./sdE/maa,'--', 'LineWidth', lw)
errorbar(sdA./sdE/sda,sesdA./sdE/sda,'--', 'LineWidth', lw)
plot(Pfgi/gst,'--', 'LineWidth', lw)
errorbar(mnT./sdE/tht,semaT./sdE/tht,'--', 'LineWidth', lw)
errorbar(mnaE./sdE/mae,semaE./sdE/mae,'--')
line([1 7], [1 1], 'Color', [.5 .5 .5], 'LineWidth', 1)
ca = gca; % ih = axes('Position',[.67 .2 .22 .4]); 
ih = axes('Position',[.4 .15 .45 .2]);  % make inset
copyobj(ca.Children, ih)
grid on, hold on, box on, axis([2.8 7.2 .965 1.025])


%%
set(gcf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, 'img/figErrorbarGAET', '-depsc', '-tiff')

