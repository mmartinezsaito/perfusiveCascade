
setPlotDefaults

% load from stats called variables

sz=4;
colcon = (exp(16*(1:-.1:0))-1)/(exp(16)-1);

figure
ha = tight_subplot(2,4,[.05],[.06 .03],[.05 .03]);
for i=1:4, h=eval(hcell{i}); hv=h.Values;
    axes(ha(i))
    [mx,my]=meshgrid(bx{i},by{i}); 
    surf(mx,my,hv'), caxis(pcsj) 
    ca=gca; ca.XLim=axi{i}(1:2); ca.YLim=axi{i}(3:4);
    if i==1, zlabel('Probability'), end
    if i==4, colorbar, end
    if i==1 || i==4, xlabel('$a_t$','Interpreter', 'latex'), else, xlabel('$k\varepsilon_t$','Interpreter', 'latex'), end
    if i==1 || i==3, ylabel('$a_{t+1}/k$','Interpreter', 'latex'), else, ylabel('$\varepsilon_{t+1}$','Interpreter', 'latex'), end
    axes(ha(4+i)), axis(axi{i}), hold on
    contourf(mx,my,hv',colcon), grid on, box on
    colormap(flipud( (exp(16*gray)-1)/(exp(16)-1) )) 
    if i==1 || i==3, ylabel('$a_{t+1}/k$','Interpreter', 'latex'), else, ylabel('$\varepsilon_{t+1}$','Interpreter', 'latex'), end  
    if i==1 || i==4, xlabel('$a_t$','Interpreter', 'latex'), else, xlabel('$k\varepsilon_t$','Interpreter', 'latex'), end
end


%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figTransMatJoint')


