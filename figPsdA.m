
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

%% closed-form expressions
xylim = [1e-3 100 1e-4 100];

% energy signal variance
mE = [sqrt(2/pi) mnaE(1:end-1)];
sde = ones(1,nl); %[1 sdE(1:end-1)];
varec = sde .* (2/3)^4;
varec = (sde .* (1-g)).^2  % because g fraction of the time, A=0 due to reset
%mtg = mnT .* (1-g).^(1/4);

Saa = @(l,f) 2 * sde(l)^2 * tr *mnT(l)^4 ./ ((2*pi*f*tr*mnT(l)^2).^2 + 1);
Saa = @(l,f) 2 * varec(l) * tr *mnT(l)^4 ./ ((2*pi*f*tr*mnT(l)^2).^2 + 1);
%Saa = @(l,f) 2 * tr * (1-g(l))^2 * (g(l)^(1/8)*mnT(l))^4 ./ ((2*pi*f*tr*(g(l)^(1/8)*mnT(l))^2).^2 + 1);
%Saa = @(l,f) 2 * tr * (1-g(l))^2 * (.9*mnT(l))^4 ./ ((2*pi*f*tr*(.9*mnT(l))^2).^2 + 1);
%Saa = @(l,f) 2 * tr * (1-g(l))^2 * mtg(l)^4 ./ ((2*pi*f*tr*mtg(l)^2).^2 + 1);
Saa = @(l,f) 2 * tr * (1-g(l))^2 * .75^2*mnT(l)^4 ./ ((2*pi*f*tr*.75*mnT(l)^2).^2 + 1);

% approximate .81 using whole integral? with E[1/T]
%fethn2numn = int(subs(fethn2,t,y)/y, y, [0.001, inf]);
%for i=1:nl, numn(i) = vpa(subs(fethn2numn,{th, n},{mnT(i), 100})), end
%numn.^-1 ./ mnT.^2 = 0.54587203

% corner freq
corfrq = 1./(2*pi*tr * 3/4 .* mnT.^2)
% freq at zero
S0 = 2 * (1-g).^2 * tr * 9/16 .* mnT.^4

x=0:1e-4:100;
figure
for  i=1:nl
    plot(x, Saa(i,x)), axis(xylim)
    hold on
end, %plot(x, sumStt(x, ts(1:end-1)))
set(gca, 'yscale', 'log', 'xscale', 'log'), grid on

%% fig04
setPlotDefaults
x=0:1e-4:100;
cocm = lines(nl);
yl = sA(:,end) - mean(sA(:,end));
xylim = [1e-3 100 1e-4 100];

fh1 = figure, fh2 = figure
for  i=1:nl
    figure(fh1)
      
    % dB
    plot(x, Saa(i,x), 'LineWidth', lw/2, 'Color', cocm(i,:))
    axis(xylim), hold on      
    
    % use squared-waved, full, data
    y = sA(:,i) - mean(sA(:,i)); 
    ns = length(y);
    winn = ns/100; % ns/8, ns/30, ns/100, ns/300, sqrt(ns)

    %[Pyy, F, Pyyc] = periodogram(y, rectwin(ns), [], fs(1), 'twosided');
    [Pyy, F, Pyyc] = pwelch(y, winn, [], [], fs(1), 'twosided');
    %[Pyy, F] = pmtm(y(1:1e6), 3, [], fs(1), 'twosided');
    %[Pyy, F, Pyyc] = pburg(y, 1, x, fs(1), 'twosided');
    
    % OR use subsampled data
    %y = sAgi{i} - mean(sAgi{i});       
    %winn = length(sAgi{i})/100;    
    %[Pyy, F, Pyyc] = pwelch(y, winn, [], [], fs(i), 'twosided', 'mean'); 
    
    ph(i) = plot(F, Pyy, ':', 'MarkerSize', msz, 'LineWidth', lw, 'Color', cocm(i,:));
    hold on
    
    % fit Lorentzian curve to PSD
    %Fft = [F(65537:end)-200; F(1:65536)]; Pyyft = [Pyy(65537:end); Pyy(1:65536)];
    %fitresult = createFit(Fft, Pyyft)
        
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
winn = size(sA,1)/100; ys = sum(sA,2) - mean(sum(sA,2)); 
[Pyy, F, Pyyc] = pwelch(ys, winn, [], [], fs(1), 'twosided');
%[Pyy, F] = pmtm(y(1:1e6), 3, [], fs(1), 'twosided');
ph(nl+1) = plot(F, Pyy, 'k:', 'MarkerSize', msz, 'LineWidth', lw);

ca = gca; 
set(ca, 'xscale', 'log', 'yscale', 'log', 'FontSize', fsz)
axis(xylim), grid on
xlabel('Frequency (Hz)'), ylabel('Power spectral density')
legend(ph, {'$a_1$','$a_2$','$a_3$','$a_4$','$a_5$','$a_6$','$a_7$', '$\sum_ia_i$'}, 'Interpreter', 'latex')




%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/fig04')


