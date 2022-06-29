
tr = 1/sif;
ts = 1 ./ (sif * [1 Pf]); % ts = 1 ./ (sif * Pf);
fs = 1 ./ ts;
w = 0.01;

% calculate timeseries at equal sampling and input rate 
% the input rate is different at each level
sAgi = {}; sTgi = {}; 
for i=1:nl
  %isgi = [true; diff(sA(:,i))~=0];% only counts changes in A
  if i==1, isgi = true(ns,1); 
  else,    isgi = sE(:,i-1)~=0; 
  end
  sAgi{i} = sA(isgi,i);
  sTgi{i} = sT(isgi,i);
end

%% fig03
setPlotDefaults
x=0:1e-4:100;
cocm = lines(nl);
yl = sE(:,end) - mean(sE(:,end));



fh1 = figure, fh2 = figure
for  i=1:nl
    figure(fh1)
      
    % use squared-waved, full, data
    winn = ns/100; % ns/8, ns/30, ns/100, ns/300, sqrt(ns)
    y = sE(:,i) - mean(sE(:,i));
    %[Pyy, F, Pyyc] = periodogram(y, rectwin(ns), [], fs(1), 'twosided');
    [Pyy, F, Pyyc] = pwelch(y, winn, [], [], fs(1), 'twosided');
    %[Pyy, F] = pmtm(y, 3, [], fs(1), 'twosided');
    %[Pyy, F, Pyyc] = pburg(y, 1, x, fs(1), 'twosided');
    
    % OR use subsampled data
    %winn = length(Y{i})/100; 
    %y = Y{i} - mean(Y{i});  
    %[Pyy, F, Pyyc] = pwelch(y, winn, [], [], fs(i), 'twosided', 'mean'); 
    
    % CPSD
    %[Pyy, F] = cpsd(y, yl, winn, [], [], fs(1), 'twosided');
    
    ph(i) = plot(F, Pyy, ':', 'MarkerSize', msz, 'LineWidth', lw, 'Color', cocm(i,:));
    %plot(F, Pyyc, ':')
    hold on

    % check the quiescent periods effect 
    %
    figure(fh2)
    subplot(nl,2,2*i-1)
    [Pyy, F] = cpsd(y, yl, winn, [], [], fs(1), 'twosided');
    plot(F, angle(Pyy), '-', 'MarkerSize', msz, 'LineWidth', alw)
    if i==1, title('CS phase (rad)'), end
    if i~=nl, set(gca, 'XTickLabel', []), end
    ylabel(['l = ' num2str(i)]), 
    set(gca, 'xscale', 'log', 'FontSize', fsz)
    axis([1e-3 100 -pi pi]), grid on, if i==nl, xlabel('Hz'), end
    subplot(nl,2,2*i)  %Cxy = (abs(Pxy).^2)./(Pxx.*Pyy)
    [Pyy, F] = mscohere(y, yl, winn, [], [], fs(1), 'twosided');
    plot(F, Pyy, '-', 'MarkerSize', msz, 'LineWidth', lw)
    if i==1, title('Coherence'), end
    if i~=nl, set(gca, 'XTickLabel', []), end
    set(gca, 'xscale', 'log', 'FontSize', fsz)
    axis([1e-3 100 0 1]), grid on, if i==nl, xlabel('Hz'), end
    %}
end
figure(fh1)
%plot(x, sumCStt(x,ts(1:end-1))) hard to calculate
winn = size(sE,1)/100; ys = sum(sE,2) - mean(sum(sE,2)); 
[Pyy, F, Pyyc] = pwelch(ys, winn, [], [], fs(1), 'twosided');
%[Pyy, F] = pmtm(y(1:1e6), 3, [], fs(1), 'twosided');
ph(nl+1) = plot(F, Pyy, 'k:', 'MarkerSize', msz, 'LineWidth', lw)

[Pyy, F, Pyyc] = pwelch(normrnd(0,1,[ns 1]), winn, [], [], fs(1), 'twosided');
%[Pyy, F] = pmtm(y(1:1e6), 3, [], fs(1), 'twosided');
ph(nl+2) = plot(F, Pyy, ':', 'MarkerSize', msz, 'LineWidth', lw, 'Color', [.3 .3 .3])

ca = gca; 
set(ca, 'xscale', 'log', 'yscale', 'log', 'FontSize', fsz)
axis([1e-3 1e2 1e-3 1e0]), grid on
xlabel('Frequency (Hz)'), ylabel('Power spectral density')
legend(ph, {'$\varepsilon_1$','$\varepsilon_2$','$\varepsilon_3$','$\varepsilon_4$','$\varepsilon_5$','$\varepsilon_6$','$\varepsilon_7$', '$\sum_i\varepsilon_i$'}, 'Interpreter', 'latex')




%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figPsdE')


