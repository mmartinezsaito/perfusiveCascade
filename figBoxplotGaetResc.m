
setPlotDefaults

xtl1 = {'$\sigma_A / \sigma_{A_*}$' '$g / g_*$' '$\langle \Theta \rangle / \theta_*$' '$\langle |\mathcal{E}|\rangle / \langle |\mathcal{E}_*| \rangle$'};
xtl2 = {'-1','-2','-3', '-4'};
xtl3 = {'3','4','5', '5.3'};
xtl = [xtl1; xtl2]';

p75 = norminv(.75); p95 = norminv(.95);
w95 = (p95-p75)/(2*p75); p75+w95*2*p75;

figure
ha = tight_subplot(2,1,[.08],[.08 .03],[.11 .03]);
for i=1:2
  if i==1
    fn1 = 'sim_sfl/rs_H0.50_ni1e5wi1e5_mw1e+00_w1e-01.mat';
    fn2 = 'sim_sfl/rs_H0.50_ni1e5wi1e5_mw1e+00_w1e-02.mat';
    fn3 = 'sim_sfl/rs_H0.50_ni1e5wi1e5_mw1e+00_w1e-03.mat';
    fn4 = 'sim_sfl/rs_H0.50_ni1e5wi1e5_mw1e+00_w1e-04.mat';
  else
    fn1 = 'sim_sfl/rs_H0.50_ni1e5_wi1e+03_mean.mat';
    fn2 = 'sim_sfl/rs_H0.50_ni1e5_wi1e+04_mean.mat';
    fn3 = 'sim_sfl/rs_H0.50_ni1e5_wi1e+05_mean.mat';
    fn4 = 'sim_sfl/rs_H0.50_ni1e5_wi2e+05_mean.mat';
  end

  load(fn1), mnace = mean(abs(ce)); sdce=sqrt(mean(ce.^2)); ns(1) = length(wmnt); 
  Y = [wsda/sdce/sda wmng/gst wmnt/sdce/tht wmne/sdce/mae]'; 
  load(fn2), mnace = mean(abs(ce)); sdce=sqrt(mean(ce.^2)); ns(2) = length(wmnt); 
  Y = [Y; [wsda/sdce/sda wmng/gst wmnt/sdce/tht wmne/sdce/mae]'];
  load(fn3), mnace = mean(abs(ce)); sdce=sqrt(mean(ce.^2)); ns(3) = length(wmnt); 
  Y = [Y; [wsda/sdce/sda wmng/gst wmnt/sdce/tht wmne/sdce/mae]'];
  load(fn4), mnace = mean(abs(ce)); sdce=sqrt(mean(ce.^2)); ns(4) = length(wmnt); 
  Y = [Y; [wsda/sdce/sda wmng/gst wmnt/sdce/tht wmne/sdce/mae]'];

  g1 = [kron(1:4,ones(1,ns(1))) kron(1:4,ones(1,ns(2))) kron(1:4,ones(1,ns(3))) kron(1:4,ones(1,ns(4)))]'; 
  g2 = [ones(4*ns(1),1); 2*ones(4*ns(2),1); 3*ones(4*ns(3),1); 4*ones(4*ns(4),1)]; 
  G1 = [repmat(xtl1(1),1,ns(1)) repmat(xtl1(2),1,ns(1)) repmat(xtl1(3),1,ns(1)) repmat(xtl1(4),1,ns(1))...
        repmat(xtl1(1),1,ns(2)) repmat(xtl1(2),1,ns(2)) repmat(xtl1(3),1,ns(2)) repmat(xtl1(4),1,ns(2))...
        repmat(xtl1(1),1,ns(3)) repmat(xtl1(2),1,ns(3)) repmat(xtl1(3),1,ns(3)) repmat(xtl1(4),1,ns(3))...
        repmat(xtl1(1),1,ns(4)) repmat(xtl1(2),1,ns(4)) repmat(xtl1(3),1,ns(4)) repmat(xtl1(4),1,ns(4))]'; 
  G2 = [repmat(xtl2(1),4*ns(1),1); repmat(xtl2(2),4*ns(2),1); repmat(xtl2(3),4*ns(3),1); repmat(xtl2(4),4*ns(4),1)]; 
  
  axes(ha(i))
  boxplot(Y,{G1 G2}, 'factorgap', 0, 'factorseparator', 1, ...
      'symbol', '', 'whisker', w95)% ,'labelverbosity', 'minor')
  line(0:17,ones(1,18), 'LineWidth',1,'LineStyle','-','Color','k')
  hold on, ymn = grpstats(Y,{g1 g2});
  plot(1:16, ymn, 'k.', 'MarkerSize', 16)
  
  set(ha(i),'YGrid','on', 'YLim',[.97 1.04], 'XTick',1:16, ... 
      'FontSize',fsz,'TickLabelInterpreter','tex')
  ylabel('Normalized activity / discharge rate')
  text([2 6.5 10 13.5], 1.035*ones(1,4), xtl1, 'Interpreter', 'latex', 'FontSize', fsz*1.2)
  if i==1
    set(ha(i),'XTickLabels',xtl2), xlabel('Learning rate w')
    % MOVE MANUALLY
    xlp = ha(i).XLabel.Position; xlp(2) = xlp(2)-20;
    xlabel('Learning rate ($\log_{10}$ w)', 'Position', xlp, 'Interpreter', 'latex')  
  else
    set(ha(i),'XTickLabels',xtl3), xlabel('Pool size')    
    % MOVE MANUALLY
    xlp = ha(i).XLabel.Position; xlp(2) = xlp(2)-20;
    xlabel('Pool size ($\log_{10}$ N)', 'Position', xlp, 'Interpreter', 'latex')  
  end
  
  obj = findobj(ha(i), 'type', 'line');    set(obj, 'Color', 'k')
  obj = findobj(ha(i), 'tag', 'Outliers'); set(obj, 'MarkerEdgeColor', 'k')
end



%%
set(gcf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, 'img/figErrorbarGaetResc', '-depsc', '-tiff')
