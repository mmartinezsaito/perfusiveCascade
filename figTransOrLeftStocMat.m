
setPlotDefaults

% load from stats called variables

sz=4; 
colcon = (exp(8*(1:-.1:0))-1)/(exp(8)-1);

figure
ha = tight_subplot(3,4,[.05],[.06 .03],[.05 .03]);
for i=1:4
    axes(ha(i))
    [mx,my]=meshgrid(bx{i},by{i}); 
    surf(mx,my,hcv{i}'), caxis([0.1 .5])  %mesh
    %colormap(flipud(gray)), colorbar
    ca=gca; ca.XLim=axi{i}(1:2); ca.YLim=axi{i}(3:4);
    if i==1, zlabel('Probability'), end
    if i==1 || i==4, xlabel('$a_t$','Interpreter', 'latex'), else, xlabel('$k\varepsilon_t$','Interpreter', 'latex'), end
    if i==1 || i==3, ylabel('$a_{t+1}/k$','Interpreter', 'latex'), else, ylabel('$\varepsilon_{t+1}$','Interpreter', 'latex'), end
    axes(ha(4+i)), axis(axi{i}), hold on
    contourf(mx,my,hcv{i}',colcon), grid on, box on
    colormap(flipud( (exp(8*gray)-1)/(exp(8)-1))) 
    %if i==1 || i==4, xlabel('$a_t$','Interpreter', 'latex'), else, xlabel('$k\varepsilon_t$','Interpreter', 'latex'), end
    if i==1 || i==3, ylabel('$a_{t+1}/k$','Interpreter', 'latex'), else, ylabel('$\varepsilon_{t+1}$','Interpreter', 'latex'), end
    axes(ha(8+i)), hold on, grid on, axis(axi{i})
    plot(bx{i},condabxpc{i}, 'k-', 'LineWidth',1.5)
    plot(bx{i},condabxpc{i}.^2, 'k:', 'LineWidth',1.5)
    plot(bx{i},condxpc{i},'k--', 'LineWidth',1.5), box on
    if i==1 || i==4, xlabel('$a_t$','Interpreter', 'latex'), else, xlabel('$k\varepsilon_t$','Interpreter', 'latex'), end
    if i==1 || i==3, ylabel('$a_{t+1}/k$','Interpreter', 'latex'), else, ylabel('$\varepsilon_{t+1}$','Interpreter', 'latex'), end
    if i==3, legend({'$\langle |x| \rangle$','$\langle x^2 \rangle$','$\langle x \rangle$'},'Location','south','Interpreter','latex'), end
    if i==4, colorbar('South'), caxis(pcsc), end
end


%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figTransMat')


