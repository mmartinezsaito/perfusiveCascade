
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


%% Fs/2=NyquistFreq;  Fs=1 and range='half' yields power in [0, Fs/2]; 
% freqres~10kHz, which is much lower than 1e7Hz for input length 1e7 (1e5Hz if Fs=200)
%  the reason is Welch's methods poorer resolution in exchange of lower noise
% figure, periodogram(sum(T,2),[],[],sif),set(gca, 'yscale', 'log', 'xscale', 'log') % check that resolution is 1e7Hz
y = sum(sT,2)-mean(sum(sT,2));
figure
subplot(221)
pwelch(y, sqrt(ns), [], [], sif, 'onesided', 'mean')%, 'ConfidenceLevel',.95)
set(gca, 'xscale', 'log')

%% pwelch_pl
% windowdefault=sqrt(ns*.9);  overlap=0.5;  Nfft=window;  Fs=sif;  conf=0.95
% range=half;  plot_type=loglog;  results=power
figure
for i=1:nl
  pwelch_pl(sT(:,i)-mnT(i), [], [], [], sif, .95, 'half', 'loglog', 'long-mean')
  hold on
end, grid on, %set(gca, 'yscale', 'log', 'xscale', 'log')
pwelch_pl(sum(sT,2)-sum(mnT), [], [], [], sif, .95)
title('T'), xlabel('Hz'), ylabel('PSD')
legend('L1','L2','L3','L4','L5','L6','L7','sum')

% If the mean is not removed from the signal there is a large spectral
%  peak at zero frequency and the sidelobes of this peak are likely to
%  swamp the rest of the spectrum.  For this reason, the default behaviour
%  is to remove the mean.  However, the matlab pwelch does not do this.

% Pxx is the distribution of power per unit frequency. For real signals,
%  pwelch returns the one-sided PSD by default; for complex signals, it
%  returns the two-sided PSD.  Note that a one-sided PSD contains the
%  total power of the input signal

% For a length-N fft, the proper normalization is
%  Y=fft(X)/sqrt(N);
%  XX=sqrt(N)*ifft(Y);

%% pwelch (MATLAB)
% plot is done with pow2db, ie 10*log10()
% windowdefault=ns/8; noverlap=0.5;  Fs=sif;

% spectra: sT(square waves), sTgi(input-induced changes)
Y = sTgi; % sT, sTgi

%the PSD is totally different for fs>2 and fs<2. why?


colmap = lines(nl+1);
figure
for i=1:nl
  % use squared-waved, full, data
  winn = ns/100; % ns/8, ns/30, ns/100, ns/300, sqrt(ns)
  y = sT(:,i) - mean(sT(:,i));
  pwelch(y, winn, [], [], fs(1), 'twosided', 'mean')
  % OR use subsampled data
  %winn = length(Y{i})/100; 
  %y = Y{i} - mean(Y{i});  
  %pwelch(y, winn, [], [], fs(i), 'twosided', 'mean')%, 'ConfidenceLevel',.95)
  
  hold on
end
winn = size(sT,1)/300;
y = sum(sT,2) - mean(sum(sT,2)); 
pwelch(y, winn, [], [], sif, 'centered', 'mean', 'ConfidenceLevel', .95)

setPlotDefaults
ca = gca; grid on
set(ca, 'xscale', 'log', 'FontSize', fsz, 'LineWidth', alw)
%axis([-100 100 1e-8 1e2])
for i=1:nl+1
    set(ca.Children(i+2), 'Color', colmap(i,:), 'LineWidth', lw, 'MarkerSize', msz); 
end
xlim([1e-3 1e2])
title('T'), xlabel('Hz'), %ylabel('PSD')
legend('L1','L2','L3','L4','L5','L6','L7','sum')

%% closed-form expressions
w ./ (ts*2*pi) % w/tr: corner frequency
(0.53*mnT).^2 .* ts(1:end-1) % peak

% var|A+E|
sdaATdivmnTzl = [0.7557,0.5664,0.5393,0.5301,0.5395,0.5253,0.5452];
sdaATdivmnT   = [0.7558,0.5690,0.4913,0.4358,0.4092,0.3919,0.3845];


xylim = [1e-3 1e2 1e-8 1e2];
Stt = @(l,f,ts) (mnT(l) * 0.53).^2 * (2*pi)^-1.5 * ts ./ (1 + (ts/w*f).^2); % unitary
Stt = @(l,f,ts) (mnT(l) * 0.53).^2 * ts ./ (1 + (ts/w*f).^2); % non-unitary
Stt = @(l,f,ts) (mnT(l) * 0.53).^2 * ts ./ (1 + (ts/w * 2*pi*f).^2); % non-unitary, ordinary frequency
% non-unitary, ordinary frequency, with numerical sdaA
Stt = @(l,f,ts) (mnT(l).*sdaATdivmnTzl(l)).^2 * ts ./ (1 + (ts/w * 2*pi*f).^2);  
% sum of all PSDs + 2*Re(sum of all pairwise CPSDs)
sumStt = @(f, ts) sum((mnT .* sdaATdivmnTzl).^2 .* ts ./ (1 + (2*pi*f' * ts/w).^2),2);  
%sumCStt = @(f, ts) sumStt(f,ts) + ...

x=0:1e-4:100;
figure
for i=1:nl
    % dB
    plot(x, Stt(i,x,ts(i)))
    axis(xylim), hold on
end, plot(x, sumStt(x, ts(1:end-1)))
set(gca, 'yscale', 'log', 'xscale', 'log'), grid on

%% fig02
setPlotDefaults
x=0:1e-4:100;
cocm = lines(nl);
yl = sT(:,end) - mean(sT(:,end));

fh1 = figure, fh2 = figure, fh3 = figure
for  i=1:nl
    figure(fh1)
    % dB
    %ezplot(@(f) Stt(i,f,ts(i)), xylim)
    plot(x, Stt(i,x,ts(i)), 'LineWidth', lw/2, 'Color', cocm(i,:))
    axis(xylim), hold on
    
    % use squared-waved, full, data
    winn = ns/100; % ns/8, ns/30, ns/100, ns/300, sqrt(ns)
    y = sT(:,i) - mean(sT(:,i));
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

    %{
    figure(fh2)
    subplot(nl,1,i)
    plot(F, angle(Pyy), '-', 'MarkerSize', msz, 'LineWidth', alw)
    set(gca, 'xscale', 'log', 'FontSize', fsz)
    axis([1e-3 100 -pi pi]), grid on
    if i==1, title('CS phase (rad)'), end
    if i~=nl, set(gca, 'XTickLabel', []), end
    ylabel(['L' num2str(i)]), xlabel('Hz')
    hold on
    figure(fh3)
    subplot(nl,1,i)
    [Pyy, F] = mscohere(y, yl, winn, [], [], fs(1), 'twosided');
    plot(F, Pyy, '-', 'MarkerSize', msz, 'LineWidth', lw)
    set(gca, 'xscale', 'log', 'FontSize', fsz)
    axis([1e-3 100 0 .3]), grid on
    if i==1, title('Magnitude-squared Coherence Estimate'), end
    if i~=nl, set(gca, 'XTickLabel', []), end
    ylabel(['L' num2str(i)]), xlabel('Hz')
    hold on
    %}
end
figure(fh1)
%plot(x, sumCStt(x,ts(1:end-1))) hard to calculate
winn = size(sT,1)/100; ys = sum(sT,2) - mean(sum(sT,2)); 
[Pyy, F, Pyyc] = pwelch(ys, winn, [], [], fs(1), 'twosided');
%[Pyy, F] = pmtm(y(1:1e6), 3, [], fs(1), 'twosided');
ph(nl+1) = plot(F, Pyy, 'k:', 'MarkerSize', msz, 'LineWidth', lw)

ca = gca; 
set(ca, 'xscale', 'log', 'yscale', 'log', 'FontSize', fsz)
axis([1e-3 100 1e-8 1e2]), grid on
xlabel('Frequency (Hz)'), ylabel('Power spectral density')
legend(ph, {'$\theta_1$','$\theta_2$','$\theta_3$','$\theta_4$','$\theta_5$','$\theta_6$','$\theta_7$', '$\sum_i\theta_i$'}, 'Interpreter', 'latex')
% make inset
ih = axes('Position',[.2 .2 .2 .2]); 
copyobj(ca.Children, ih)
set(ih, 'xscale', 'log', 'yscale', 'log')
axis([1 1.5 5e-5 2e-4]), box on, grid on


%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figPsdT')



%{
  
% PSD scaling exponents I will need help with the points indicated by yellow-highlighted dashes.

% T    b=2
figure
for i=1:nl
  pwelch_pl(T(:,i)-mnT(i), [], [], [], sif, .95, 'half', 'loglog', 'long-mean')
  hold on
end, grid on, %set(gca, 'yscale', 'log', 'xscale', 'log')
pwelch_pl(sum(T,2)-sum(mnT), [], [], [], sif, .95)
title('T'), xlabel('Hz'), ylabel('PSD')
legend('L1','L2','L3','L4','L5','L6','L7','sum')
preprint

% A   b=1.5 % pwelch_pl
figure
subplot(121)
for i=1:nl, pwelch_pl(A(:,i), [], [], [], sif), hold on, end, grid on
pwelch_pl(sum(A,2), [], [], [], sif), set(gca, 'yscale', 'log', 'xscale', 'log')
legend('L1','L2','L3','L4','L5','L6','L7','sum'), title('A')
subplot(122)
for i=1:nl, pwelch_pl(Abf(:,i), [], [], [], sif), hold on, end, grid on
pwelch_pl(sum(Abf,2), [], [], [], sif), set(gca, 'yscale', 'log', 'xscale', 'log')
legend('L1','L2','L3','L4','L5','L6','L7','sum'), title('Abf')
preprint

% E    b=0
figure
subplot(121), for i=1:nl, pwelch_pl(E(:,i), [], [], [], sif), hold on, end, grid on
pwelch_pl(sum(E,2), [], [], [], sif), set(gca, 'yscale', 'log', 'xscale', 'log')
legend('L1','L2','L3','L4','L5','L6','L7','sum'), title('E')
subplot(122)
for i=1:nl, pwelch_pl(Ebf(:,i), [], [], [], sif), hold on, end, grid on
pwelch_pl(abs(sum(Ebf,2)), [], [], [], sif), set(gca, 'yscale', 'log', 'xscale', 'log')
title('Ebf'), legend('L1','L2','L3','L4','L5','L6','L7','sum')
preprint

% td+bu
figure
subplot(121), for i=1:nl, pwelch_pl(E(:,i)+A(:,i), [], [], [], sif), hold on, end
pwelch_pl(sum(E+A,2), [], [], [], sif), set(gca, 'yscale', 'log', 'xscale', 'log')
title('L+E'), legend('L1','L2','L3','L4','L5','L6','L7','sum'), grid on
subplot(122)
for i=1:nl, pwelch_pl(Ebf(:,i)+Abf(:,i), [], [], [], sif), hold on, end
pwelch_pl(abs(sum(Ebf+Abf,2)), [], [], [], sif), set(gca, 'yscale', 'log', 'xscale', 'log'),
title('Lbf+Ebf'), legend('L1','L2','L3','L4','L5','L6','L7','sum'), grid on
 preprint

%}
