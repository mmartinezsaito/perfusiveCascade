 
setPlotDefaults %fsz=10;
lw = 0.3;

sif = 200;
tr = 1/sif;
ts = 1 ./ (sif * [1 Pf]); % ts = 1 ./ (sif * Pf);
fs = 1 ./ ts;
w = 0.01;

x=0:1e-4:100;
xylimT = [1e-3 1e2 1e-8 1e2];
xylimE = [1e-3 1e2 1e-3 1e0];
xylimA = [1e-3 1e2 1e-4 1e2]; 
xylim = xylimA

figure
%ah = tight_subplot(4,4,[.02 .03],[.06 .04],[.06 .02]);
ah = tight_subplot(4,4,[.02 .02],[.06 .04],[.05 .03]);
for i=1:2:nl  %1:nl
  %yi = sA(:,i) - mean(sA(:,i)); 
  yi = sE(:,i) - mean(sE(:,i)); 
  ns = length(yi);
  winn = ns/100; % ns/8, ns/30, ns/100, ns/300, sqrt(ns)
    
  for j=1:2:nl  %1:nl
    ca = ah(2*(i-1)+(j+1)/2);
    axes(ca) %nl*(i-1)+j %subplot(nl,nl,nl*(i-1)+j)
    yj = sA(:,j) - mean(sA(:,j));
    %yj = sE(:,j) - mean(sE(:,j));
    %{
    if i==j % single variable CPSD 
      [Pyy, F, Pyyc] = pwelch(yi, winn, [], [], fs(1), 'twosided');
      ph = plot(F, Pyy, '-', 'LineWidth', lw, 'Color', 'r');
      hold on
      [Pyy, F, Pyyc] = pwelch(yj, winn, [], [], fs(1), 'twosided');
      ph = plot(F, Pyy, '-', 'LineWidth', lw, 'Color', 'k');
      axis(xylim), set(gca, 'yscale', 'log')
    elseif i<j 
      [Pyy, F] = cpsd(yi, yj, winn, [], [], fs(1), 'twosided');
      plot(F, angle(Pyy), '-', 'MarkerSize', msz, 'LineWidth', alw)
      axis([1e-3 1e2 -pi pi])
    elseif i>j 
      [Pyy, F] = mscohere(yi, yj, winn, [], [], fs(1), 'twosided');
      plot(F, Pyy, '-', 'MarkerSize', msz, 'LineWidth', lw)      
      axis([1e-3 1e2 0 1])
    end
    if ~ismember(j,[1 i i+2]), set(gca, 'YTickLabel', []), end    
    %}
    % multiple variable CPSD
    yyaxis right
    [Pyy, F] = cpsd(yi, yj, winn, [], [], fs(1), 'twosided');
    plot(F, angle(Pyy), '-', 'MarkerSize', msz, 'LineWidth', lw, 'Color', [0 1 0 .7])
    if j~=nl, set(gca, 'YTickLabel', []), end
    axis([1e-3 1e2 -pi pi]), grid on
    yyaxis left 
    [Pyy, F] = mscohere(yi, yj, winn, [], [], fs(1), 'twosided');
    plot(F, Pyy, '-', 'MarkerSize', msz, 'LineWidth', lw, 'Color', [0 0 1 1])   
    if j~=1, set(gca, 'YTickLabel', []), end
    axis([1e-3 1e2 0 1])
    %}
    set(gca, 'xscale', 'log', 'XTickLabelMode', 'auto', 'FontSize', fsz, 'LineWidth', alw)     
    %if i==1, title(['\rm l = ' num2str(j)]), end
    if i==1, title(['$a_' num2str(j) '$'],'Interpreter','latex'), end
    %if j==1, ylabel(['l = ' num2str(i)]), end
    if j==1, ylabel(['$\varepsilon_' num2str(i) '$'],'Interpreter','latex', 'Color', 'k'), end
    if i==nl, xlabel('Frequency (Hz)')
    else,     set(gca, 'XTickLabel', []), end        
    grid on, box on
  end
end

%% Analytic signal through Hilbert transform

figure
for i = 1:nl
  %yi = sE(:,i) - mean(sE(:,i)); 
  %yi = sA(:,i) - mean(sA(:,i));
  yi = sE(:,i) - mean(sE(:,i)) + sA(:,i) - mean(sA(:,i));
  h = hilbert(yi);
  env = abs(h);
  [Pyy, F] = pwelch(env, winn, [], [], fs(1), 'twosided');
  ph = plot(F, Pyy, '-', 'LineWidth', lw);
  axis(xylim), set(gca, 'xscale', 'log', 'yscale', 'log'), grid on, hold on
end    
      
%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figEpsactcoh')  %Epscoh, Actcoh, Thrcoh, Epsactcoh

