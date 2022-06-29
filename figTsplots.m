tr = 1/sif;
ts = 1 ./ (sif * [1 Pf]);
fs = 1 ./ ts;
w = 0.01;
setPlotDefaults

%% figT
figure
plot(sT), axis([5e5 5e5+5e3/tr 0 15])

%% figAE
lwl = .1; lwr = 1;
figure
[ah, pos] = tight_subplot(3,2,[.04 .04],[.06 .03],[.06 .03]);
for i=1:2
  axes(ah(0+i)) %subplot(3,2,0+i)   
  if i==1, ylabel('Pulse activity'), lw=lwl; else, ylabel(''), lw=lwr; set(gca, 'YTickLabel', []), end
  ph=plot(sE, 'LineWidth', lw);
  if i==2, legend(ph, {'$\varepsilon_1$','$\varepsilon_2$','$\varepsilon_3$','$\varepsilon_4$','$\varepsilon_5$','$\varepsilon_6$','$\varepsilon_7$'}, 'Interpreter', 'latex','orientation','horizontal','FontSize',fsz,'LineWidth',alw), end
  xlabel(''), set(gca, 'XTickLabel', [], 'FontSize', fsz, 'LineWidth', alw)
  axis([5e5 5e5+10/tr/10^(i-1) -15 15]), grid on  
  axes(ah(2+i)) %subplot(3,2,2+i)  
  ph=plot(sA, 'LineWidth', lw);
  if i==1, ylabel('Subthreshold activity'), else, ylabel(''), set(gca, 'YTickLabel', []), end
  if i==2, legend(ph, {'$a_1$','$a_2$','$a_3$','$a_4$','$a_5$','$a_6$','$a_7$'}, 'Interpreter', 'latex','orientation','horizontal','FontSize',fsz,'LineWidth',alw), end
  xlabel(''), set(gca, 'XTickLabel', [], 'FontSize', fsz, 'LineWidth', alw)
  axis([5e5 5e5+10/tr/10^(i-1) -10 10]), grid on 
  axes(ah(4+i)) %subplot(3,2,4+i)  
  phE=plot(sum(sE,2),'LineWidth', lw, 'Color', 'r'); hold on
  phA=plot(sum(sA,2),'LineWidth', lw, 'Color', 'k'); 
  phEA=plot(sum(sA+sE,2),'LineWidth', lw, 'Color', 'g');
  if i==1, ylabel('Pulse + subthreshold activity'), else, ylabel(''), set(gca, 'YTickLabel', []), end
  if i==2, legend([phE;phA;phEA], {'$\sum_i\varepsilon_i$','$\sum_ia_i$','$\sum_i\varepsilon_i+a_i$'}, 'Interpreter','latex','FontSize',fsz,'LineWidth',alw), end
  xlabel('Time (s)'), set(gca, 'FontSize', fsz, 'LineWidth', alw)
  if i==1, set(gca, 'XTickLabel', {'0' '2.5' '5' '7.5' '10'}) %{'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'})
  else, set(gca, 'XTickLabel', {'0' '0.25' '0.5' '0.75' '1'}),end %{'0' '0.1' '0.2' '0.3' '0.4' '0.5' '0.6' '0.7' '0.8' '0.9' '1'}), end
  axis([5e5 5e5+10/tr/10^(i-1) -20 20]), grid on
end

%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figAEts')   %  figabsAEts


%% figabsAplusE
figure
ah = tight_subplot(2,1,[.04 .04],[.06 .03],[.1 .03]);
for i=1:2
  axes(ah(i))%, subplot(2,1,i)
  plot(sum(abs(sE),2), 'Color', 'r'), hold on
  plot(sum(abs(sA),2), 'Color', 'k')
  plot(sum(abs(sA)+abs(sE),2), 'Color', 'g')
  ylabel('Energy'), xlabel(''), grid on
  legend({'$\sum_i|\varepsilon_i|$','$\sum_i|a_i|$','$\sum_i|\varepsilon_i|+|a_i|$'}, 'Interpreter', 'latex','FontSize',fsz,'LineWidth',alw)
  axis([4e5 4e5+10/tr/10^(i-1) 0 50])
  set(gca, 'FontSize', fsz, 'LineWidth', alw)
  if i==1, set(gca, 'XTickLabel', {'0' '2.5' '5' '7.5' '10'})
  else, set(gca, 'XTickLabel', {'0' '0.25' '0.5' '0.75' '1'}), end
end
xlabel('Time (s)')

%% figSqAplusE
figure
ah = tight_subplot(2,1,[.04 .04],[.06 .03],[.1 .03]);
for i=1:2
  axes(ah(i))%, subplot(2,1,i)
  plot(sum(sE.^2,2), 'Color', 'r'), hold on
  plot(sum(sA.^2,2), 'Color', 'k')
  plot(sum(sA.^2+sE.^2,2), 'Color', 'g')
  ylabel('Energy (squared activities)'), xlabel(''), grid on
  legend({'$\sum_i \varepsilon_i^2$','$\sum_i a_i^2$','$\sum_i \varepsilon_i^2+ a_i^2$'}, 'Interpreter', 'latex','FontSize',fsz,'LineWidth',alw)
  axis([4e5 4e5+10/tr/10^(i-1) 0 1e3])
  set(gca, 'FontSize', fsz, 'LineWidth', alw)
  if i==1, set(gca, 'XTickLabel', {'0' '2.5' '5' '7.5' '10'})
  else, set(gca, 'XTickLabel', {'0' '0.25' '0.5' '0.75' '1'}), end
end
xlabel('Time (s)')

%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figsqAEts')  % figSqAEpdf


