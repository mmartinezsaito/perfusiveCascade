setPlotDefaults

% ACT=A; A=L; L=V;
k = 1.5; s0 = 1.5;  % k=1.25
%ss = sprintf('img/fig01_input%s_nl%02u_ni%u_w%.2f_isrand%u_eam%.1f_ueam%.1f', inputpdf, nl, ni, w, iseamrand, eam, ueam);
ni = max(size(A));
gv = kron(1:nl, ones(1, 0.9*ni))';
sA = A(ni/10+1:end,:); %sAbf = Abf(ni/10+1:end,:);
sT = T(ni/10+1:end,:); rT = reshape(sT, 0.9*ni*nl, 1);
sE = E(ni/10+1:end,:); 

Pa0 = mean(sA==0, 1);
Pf = mean(sE~=0, 1); 
mdT = median(sT); mxT = max(sT); mnT = mean(sT);
sr = mean(sT(:,2:end)) ./ mean(sT(:,1:end-1));

rT = reshape(sT, 0.9*ni*nl, 1);
gamdT = fitdist(rT,'Gamma', 'By', gv); 
gevdT = fitdist(rT,'GeneralizedExtremeValue', 'By', gv); 

absA = {}; for i=1:nl, absA{i} = abs(sA(sA(:,i)~=0,i)); end
absE = {}; for i=1:nl, absE{i} = abs(sE(sE(:,i)~=0,i)); end
for i=1:nl, gpdE{i} = fitdist(absE{i}(absE{i}>mnT(i)) - mnT(i), 'GeneralizedPareto'); end
for i=1:nl, gamdE{i} = fitdist(absE{i}, 'Gamma'); end
tabsEs = sE(sE>0);
for i=1:nl, gpsr(i) = gpstat(gpdE{i}.k, gpdE{i}.sigma, gpdE{i}.theta); end
mnaE = [sqrt(2/pi) cellfun(@mean, absE)];

sgthL = cellfun(@(x) x.sigma, gpdE) ./ mnT;
xiL = cellfun(@(x) x.k, gpdE);

%% basic plots

nrows = 3;
ncols = nl;
maxabsL = ceil(max(abs(sA(:))));   
maxabsT = ceil(max(abs(sT(:))));  
maxabsE = ceil(max(abs(sE(:))));      

figure
for i = 1:nrows*ncols
    switch ncols
        case 3, evr = 3;
        case 4, evr = 2;
        case 7, evr = 1;
    end
    
    subplot(nrows, ncols, i)
    if i <= ncols   % absA
        j = (i-(1))*evr+1;
        h = histogram(absA{j}, maxabsL*2);
        title(['P(a_' num2str(j) '=0)=' sprintf('%.3f',Pa0(j))])
        
    elseif i > ncols && i <= ncols*2   % sT
        j = (i-(ncols+1))*evr+1;
        h = histogram(sT(:,j), maxabsT*2); hold on
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * gampdf(h.BinEdges, gamdT{j}.a, gamdT{j}.b))
        title({sprintf('\\theta_%d~%.3f , mnaE=%.3f', j, mnT(j), mnaE(j)), ...
            sprintf('sh=%.2f, sc=%.2f',gamdT{j}.a, gamdT{j}.b)})
        
    elseif i > 2*ncols && i <= ncols*3   % absE
        j = (i-(2*ncols+1))*evr+1;
        h = histogram(absE{j}, maxabsE*2); hold on
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gp', h.BinEdges, gpdE{j}.k, gpdE{j}.sigma, mdT(j)), '-')
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gamma', h.BinEdges, gamdE{j}.a, gamdE{j}.b), '-')
        title({sprintf('P(e_%d=0) = %.3f', j, Pf(j)), ...
            sprintf('\\xi=%.2f, \\sigma=%.2f', gpdE{j}.k, gpdE{j}.sigma)})
            %sprintf('sh=%.2f, sc=%.2f', gamdE{j}.a, gamdE{j}.b)})
        
        if j>1, xlabel({sprintf('\\theta/\\theta_{-1}=%.3f, %.3f', sr(j-1), 1+gpdE{j-1}.sigma/(1-gpdE{j-1}.k) / mnT(j-1)), ...
                sprintf('\\sigma/\\sigma_{-1}=%.3f', gpdE{j}.sigma / gpdE{j-1}.sigma), ...
                sprintf('\\sigma/\\theta=%.2f', gpdE{j}.sigma / mnT(j))}), 
        else,   xlabel(sprintf('\\sigma_1/\\theta_1=%.2f', gpdE{j}.sigma / mnT(j)))
        end
    end
    hdl = length(h.Data); hbw = h.BinWidth;
    ylim([0, hdl * hbw/k^j * 4]), if i > ncols && i <= ncols*2, ylim([0, hdl * hbw/k^j * 16]), end
    my = get(gca, 'YLim'); line([mnT; mnT], [zeros(1,nl); my(2)*ones(1,nl)], 'LineStyle', '--')
    xlim([0 s0 * k^j])     %xlim([0 2 * 2^j])
end


%%
% fig01 v1

figure
for i = 1:nrows*ncols
    switch ncols
        case 3, evr = 3;
        case 4, evr = 2;
        case 7, evr = 1;
    end
    
    subplot(nrows, ncols, i)
    if i <= ncols   % absL
        j = (i-(1))*evr+1;
        h = histogram(absA{j}, maxabsL*2);
            
    elseif i > ncols && i <= ncols*2   % sT
        j = (i-(ncols+1))*evr+1;
        h = histogram(sT(:,j), maxabsT*2); hold on
        plot(h.BinEdges, gampdf(h.BinEdges, gamdT{j}.a, gamdT{j}.b))

    elseif i > 2*ncols && i <= ncols*3   % absE
        j = (i-(2*ncols+1))*evr+1;
        h = histogram(absE{j}, maxabsE*2); hold on
        plot(h.BinEdges, pdf('gp', h.BinEdges, gpdE{j}.k, gpdE{j}.sigma, mdT(j)), '-')
       
    end
    %hbc = conv(h.BinEdges, [.5 .5], 'valid');
    %patch([hbc fliplr(hbc)], [h.Values+.1  fliplr(h.Values)-.1], .8*ones(1,3)); 
 
    set(h, 'Normalization','pdf', 'DisplayStyle','stairs', ...
        'LineWidth',2, 'LineStyle', ':')
    hdl = length(h.Data); hbw = h.BinWidth;
    ylim([0, 3.5*k^-j])
    my = get(gca, 'YLim'); line([mnT; mnT], [zeros(1,nl); my(2)*ones(1,nl)], 'Color', .8*ones(1,3))
    xlim([0 s0 * k^j])     %xlim([0 2 * 2^j])
end


%%
% fig01 v2

figure
sp = [1:2:8 2:2:8];
for i = 0:7
    smplot(4, 2, sp(i+1), 'axis', 'on')
    
    box on, ylim([0, 2*k^2 * k^-i]), xlim([0 s0 * k^i])   
    my = get(gca, 'YLim'); 
    line([mnT; mnT], [zeros(1,nl); my(2)*ones(1,nl)], 'LineWidth', 1, 'Color', .8*ones(1,3)) 
    %set(gca, 'FontSize', fsz, 'FontName','LM Roman 10', 'LineWidth', alw)
    set(gca, 'FontSize', fsz, 'LineWidth', alw)
    hold on
    
    if i==0
        plot(0:.1:3, exp(-(0:.1:3).^2/2) * sqrt(2/pi)')
        ylim([0, 3/2]), xlim([0 3]), title('Exogenous input')
        continue
    end
   
    [hv, hbe] = histcounts(absA{i},  'Normalization', 'pdf');
    hbc = conv(hbe, [.5 .5], 'valid');
    plot(hbc, hv, 'ko:', 'LineWidth', lw, 'MarkerSize', msz, 'LineWidth', lw);

    [hv, hbe] = histcounts(sT(:,i),  'Normalization', 'pdf');
    hbc = conv(hbe, [.5 .5], 'valid');
    plot(hbc, hv, 'go ', 'LineWidth', lw, 'MarkerSize', msz);
    plot(hbc, gampdf(hbc, gamdT{i}.a, gamdT{i}.b), 'g-', 'LineWidth', lw)
   
    [hv, hbe] = histcounts(absE{i},  'Normalization', 'pdf');
    hbc = conv(hbe, [.5 .5], 'valid');
    plot(hbc, hv, 'ro ', 'LineWidth', lw, 'MarkerSize', msz);
    plot(hbc, pdf('gp', hbc, gpdE{i}.k, gpdE{i}.sigma, mnT(i)), 'r-', 'LineWidth', lw)
   
    %patch([hbc fliplr(hbc)], [h.Values+.1  fliplr(h.Values)-.1], .8*ones(1,3)); 
 
    hdl = length(hv); hbw = mean(diff(hbe));
    switch i
        case 1,    title('l = 1 (sensory)')
        case 7,    title('l = 7 (top)')
        otherwise, title(sprintf('l = %u',i))
    end
end

%%

cf = gcf; ca = gca; 
% PaperSize, PaperPosition00, PaperType are for paged formats such as pdf or ps
% set(cf, 'PaperUnits','centimeters', 'PaperSize', [8.5 11],...
%    'PaperPositionMode','auto', 'PaperPosition',[0 0 8.5 11])
%  print(gcf, '-dpdf', 'img/fig01', '-bestfit')
%set(cf, 'PaperUnits','centimeters', 'PaperPositionMode','manual', 'PaperPosition',[0 0 8.5 11])
set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figAETpdf')



%{
xtl{1} = {'L2' 'L4' 'L6' 'L8'};
xtl{2} = {'L2' 'L4' 'L6' 'L8'};
ytl{1} = {'0' '0.2' '0.4' '0.6' '0.8' '1'}; 
ytl{2} = {'0.4' '' '0.6' '' '0.8' '' '1'};
fgid  = char(97:103);


for i = 5%1:nrows * ncols
    fstr = sprintf('fig02%c', fgid(i))
    fh = openfig(fstr);
    figure(fh)
    set(fh, 'Units', 'inches', 'Position', [1 1 width height]);
    if mod(i, ncols)~=1, ylabel(''), end
    if i <= ncols,       xlabel(''), end
    title('')
    fc = fh.Children;
    ca = gca;
    set(ca, 'FontSize', fsz, 'LineWidth', alw, ...
        'XTickLabel', xtl{(mod(i,ncols)==0)+1}, 'YTickLabel', ytl{(mod(i,ncols)==0)+1})
    cc = ca.Children;
    for j = 1:length(ca.Children)
        set(cc(j), 'LineWidth', lw)
        if strcmp(get(cc(j),'Type'), 'errorbar')
            set(cc(j), 'CapSize', csz)
            if i==5, set(cc(1), 'LineStyle', '--'), set(cc(3), 'LineStyle', '-'), end
        end
        if strcmp(get(cc(j),'Type'), 'line'), set(cc(j), 'LineWidth', lw), end
       
    end
    
    
    % Printing
    set(fh, 'PaperUnits', 'inches', 'PaperPosition', [0 0 width height])
    print(fh, '-depsc', fstr, '-r300')    % EPS
    print(fh, '-dtiff', fstr, '-r300')             % TIFF
end

cf = gcf; set(cf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
print(gcf, '-djpeg', [ss '_fig' sprintf('%02u', cf.Number) '.jpeg'], '-r200')
%}
