setPlotDefaults

%load sim_sfl/rs_ni1e5wi1e5_mean.mat
%load sim_sfl/rs_H0.50_ni1e5wi1e5_w1e-03.mat
load sim_sfl/rs_H0.50_ni1e5wi1e5_mw1e+00_w1e-02.mat

cg = mean(wmng), cg2=mean(ca==0)
mnace = mean(abs(ce)) % rescaled to 1
sdace = std(abs(ce))
sdce = sqrt(mean(ce.^2))  % sqrt(mnace^2+sdace^2); std(ce)
mnaca = mean(abs(ca)), sdaca = std(abs(ca)), sdca = mean(ca.^2)^.5 % sdca=std(ca)
mean(abs(ca(ca~=0)))/sdce
mnacae = mean(abs(cae)), mdcae = median(abs(cae))
cket = mnace/mnacae 
ckte = mnacae, ck = mnace      % th mean, resca mean

mnacae / sdce % 
sdca / sdce % std(ca); sqrt(mean(ca.^2))


%% figSFL
lw = 1.5; msz = 5;
ca(ca==0) = [];

figure, box on, grid on, axis([0 3 0 3]), hold on
my = get(gca, 'YLim'); 
set(gca, 'FontSize', fsz, 'FontName','LM Roman 10', 'LineWidth', alw)

%line([1;1]*[mnaca mnace mnacae sdce], [zeros(1,4); my(2)*ones(1,4)], 'LineWidth', lw, 'Color', [.8;.2;.2;.2]*ones(1,3)) 
ph(1)=line([1;1]*sdca, [0; my(2)], 'LineWidth', 1, 'Color', 'k','LineStyle','-.') 
ph(2)=line([1;1]*mnace, [0; my(2)], 'LineWidth', 1, 'Color', 'r','LineStyle','--')
ph(3)=line([1;1]*mnacae, [0; my(2)], 'LineWidth', 1, 'Color', 'g','LineStyle','--') 
ph(4)=line([1;1]*sdce, [0; my(2)], 'LineWidth', 1, 'Color', 'r','LineStyle','-.') 

[hv, hbe] = histcounts(abs(ca),  'Normalization', 'pdf'); hbc = conv(hbe, [.5 .5], 'valid');
ph(5)=plot(hbc, hv, 'ko-', 'LineWidth', lw, 'MarkerSize', msz, 'LineWidth', lw);
[hv, hbe] = histcounts(ct,  'Normalization', 'pdf'); hbc = conv(hbe, [.5 .5], 'valid');
ph(6)=plot(hbc, hv, 'go-', 'LineWidth', lw, 'MarkerSize', msz, 'LineWidth', 1);
[hv, hbe] = histcounts(abs(ce),  'Normalization', 'pdf'); hbc = conv(hbe, [.5 .5], 'valid');
ph(7)=plot(hbc, hv, 'ro-', 'LineWidth', lw, 'MarkerSize', msz, 'LineWidth', lw);
[hv, hbe] = histcounts(abs(cae),  'Normalization', 'pdf'); hbc = conv(hbe, [.5 .5], 'valid');
ph(8)=plot(hbc, hv, 'bo-', 'LineWidth', lw, 'MarkerSize', msz, 'LineWidth', 1);
%{
[hv, hbe] = histcounts(wmnt,  'Normalization', 'pdf'); hbc = conv(hbe, [.5 .5], 'valid');
plot(hbc, hv, 'g--', 'LineWidth', lw, 'MarkerSize', msz);
[hv, hbe] = histcounts(wmng,  'Normalization', 'pdf'); hbc = conv(hbe, [.5 .5], 'valid');
plot(hbc, hv, 'co--', 'LineWidth', lw, 'MarkerSize', msz, 'LineWidth', lw);
%}
legend(ph, {'$\sigma_{A^+}$','$\langle|\tilde{A}|\rangle$','$\langle|\mathcal{E}|\rangle$',...
    '$\sigma_\mathcal{E}$','$|A^+|$','$\Theta$','$|\mathcal{E}|$','$|\tilde{A}|$'}, ...
    'Interpreter', 'latex')
ylabel('Probability density'), xlabel('Rescaled activity')
set(gca,'FontName','Arial')

%%
set(gcf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, 'img/figSFL', '-depsc', '-tiff')




