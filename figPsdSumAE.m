
tr = 1/sif;
ts = 1 ./ (sif * [1 Pf]); % ts = 1 ./ (sif * Pf);
fs = 1 ./ ts;
w = 0.01;

%% figPsdsumAE
setPlotDefaults
x=0:1e-4:100;
%xylim = [1e-3 100 1e-3 100];
xylim = [1e-3 100 5e-2 5e3]; 

figure
winn = size(sA,1)/100; 
%yas = sum(abs(sA),2) - mean(sum(abs(sA),2)); 
yas = sum(sA.^2,2) - mean(sum(sA.^2,2)); 
[Pyy, F, Pyyc] = pwelch(yas, winn, [], [], fs(1), 'twosided');
ph(1) = plot(F, Pyy, 'k:', 'MarkerSize', msz, 'LineWidth', lw);
hold on

winn = size(sE,1)/100; 
%yes = sum(abs(sE),2) - mean(sum(abs(sE),2)); 
yes = sum(sE.^2,2) - mean(sum(sE.^2,2)); 
[Pyy, F, Pyyc] = pwelch(yes, winn, [], [], fs(1), 'twosided');
ph(2) = plot(F, Pyy, 'r:', 'MarkerSize', msz, 'LineWidth', lw);

winn = size(sE,1)/100; ys = yas + yes; 
[Pyy, F, Pyyc] = pwelch(ys, winn, [], [], fs(1), 'twosided');
ph(3) = plot(F, Pyy, 'g:', 'MarkerSize', msz, 'LineWidth', lw);

ca = gca; 
set(ca, 'xscale', 'log', 'yscale', 'log', 'FontSize', fsz)
axis(xylim), grid on
xlabel('Frequency (Hz)'), ylabel('Power spectral density')
%legend(ph, {'$\sum_i |a_i|$','$\sum_i |\varepsilon_i|$', '$\sum_i |a_i|+|\varepsilon_i|$'}, 'Interpreter', 'latex')
legend(ph, {'$\sum_i a_i^2$','$\sum_i \varepsilon_i^2$', '$\sum_i a_i^2+\varepsilon_i^2$'}, 'Interpreter', 'latex')


%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figPsdsumsqAE')





