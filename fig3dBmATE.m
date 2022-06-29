
setPlotDefaults
fsz = 14

%load sim_sfl/rs_H0.50_ni1e5_wi1e+05_mean.mat
load sim_sfl/rs_H0.50_ni1e6_wi2e+05_mean_v2.mat
%load sim_sfl/rs_H0.50_ni1e5wi1e5_mw1e+00_w1e-02.mat

figure
plot3(wsda/sdce,wmnt/sdce,wmne/sdce, 'k'), grid on, hold on, axis(axn)
plot3(axn(2)*ones(1,ni),wmnt/sdce,wmne/sdce, 'Color',[.5 .5 .5])
plot3(wsda/sdce,axn(4)*ones(1,ni),wmne/sdce, 'Color',[.5 .5 .5])
plot3(wsda/sdce,wmnt/sdce,axn(5)*ones(1,ni), 'Color',[.5 .5 .5])
plot3(sda,tht,mae,'ko','markersize',radius*1e3), plot3(sda,tht,mae,'k*')
plot3(axn(2),tht,mae,'k*'), plot3(sda,axn(4),mae,'k*'),plot3(sda,tht,axn(5),'k*')

xlabel('$\sigma_A / \sigma_{A_*}$', 'Interpreter', 'latex','FontSize',fsz)  
ylabel('$\langle \Theta \rangle / \theta_*$', 'Interpreter', 'latex', 'FontSize',fsz)  
zlabel('$\langle |\mathcal{E}|\rangle / \langle |\mathcal{E}_*|\rangle$', 'Interpreter', 'latex', 'Fontsize',fsz)  


%%

set(gcf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, 'img/fig3dBmATE', '-depsc', '-tiff')
