nr = ceil(sqrt(nl));

%% A
sL = A(ni/10+1:end,:); sLbf = Abf(ni/10+1:end,:);
Pa0 = mean(sL==0, 1);
absL = {}; for i=1:nl, absL{i} = abs(sL(sL(:,i)~=0,i)); end
%for i=1:nl, gevd{i} = fitdist(absL{i}, 'GeneralizedExtremeValue'); end
%for i=1:nl, gamd{i} = fitdist(absL{i}, 'Gamma'); end
tabsL = sL(sL~=0);
figure, nr=ceil(sqrt(nl+2));
for i=1:nl+2
    subplot(nr,nr,i)
    if     i == nl+1, histogram(abs(sum(tabsL,2))), title(['abssumL, P(a_{all}=0)=' num2str(mean(Pa0))])
    elseif i == nl+2, histogram(sum(abs(tabsL),2)), title('sumabsL')
    else,  histogram(absL{i}) 
        title(['P(a_' num2str(i) '),  P(a_' num2str(i) '=0)=' num2str(Pa0(i))])
    end, %set(gca, 'yscale', 'log')    
end

%% T
gv = kron(1:nl, ones(1, 0.9*ni))';
sT = T(ni/10+1:end,:); 
rT = reshape(sT, 0.9*ni*nl, 1);
gamd = fitdist(rT,'Gamma', 'By', gv); 
gevd = fitdist(rT,'GeneralizedExtremeValue', 'By', gv); 
figure
for i=1:nl+2
    subplot(nr,nr,i)
    if     i == nl+1, histogram(sum(T,2)), title('sumT')
    elseif i == nl+2, histogram(sum(abs(T),2)), title('sumabsT')
    else,         h = histogram(T(:,i));
        hold on, title(['\theta_' num2str(i) '=' num2str(mdT(i)) '; shape=' num2str(gamd{i}.a) ', scale=' num2str(gamd{i}.b)])
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * gampdf(h.BinEdges, gamd{i}.a, gamd{i}.b))
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gev', h.BinEdges, gevd{i}.k, gevd{i}.sigma, gevd{i}.mu), '--')
    end, %set(gca, 'yscale', 'log')        
end

%% E

% levelwise plots
%sE = E(ni/10+1:end,:); Pf = mean(sE~=0, 1); 
absE = {}; for i=1:nl, absE{i} = abs(sE(sE(:,i)~=0,i)); end
%for i=1:nl, gevd{i} = fitdist(absE{i}, 'GeneralizedExtremeValue'); end
%for i=1:nl, gamd{i} = fitdist(absE{i}, 'Gamma'); end
for i=1:nl
    absE2fit{i} = absE{i}(absE{i}>mnT(i));
    gpdd{i} = fitdist(absE2fit{i}, 'GeneralizedPareto', 'theta', mnT(i)); 
end
tabsE = sE(sE~=0);
figure
for i=1:nl+2
    subplot(nr,nr,i)
    if     i == nl+1, histogram(abs(sum(tabsE,2))), title('abssumE')
    elseif i == nl+2, histogram(sum(abs(tabsE),2)), title('sumabsE')
    else,         h = histogram(absE{i});
        hold on, title(['P(e_' num2str(i) '),  Pf=' num2str(Pf(i))])
        %plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gev', h.BinEdges, gevd{i}.k, gevd{i}.sigma, gevd{i}.mu), '--')
        %plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gam', h.BinEdges, gamd{i}.a, gamd{i}.b), '-')
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gp', h.BinEdges, gpdd{i}.k, gpdd{i}.sigma, gpdd{i}.theta))
    end, %set(gca, 'yscale', 'log')    
end

% E is not amount of units, but amount of surprisal relayed
dv1 = diff(L==1);
if dv1(find(dv1~=0, 1, 'last'))==1,      dv1 = [1 dv1 -1]; 
elseif dv1(find(dv1~=0, 1, 'last'))==-1, dv1 = [1 dv1]; 
end
avlen = find(dv1 == -1) - find(dv1 == 1);

% activities increase as ~ k
% probabilities decrease as ~ g ~ k^-2
% so the exponent is  ~ -2 for local fields
%    sumaE ~ k^2, so  ~ -1 for global fields
% local fields
mnaE
diff(log(mnaE))
% global summed field
sumaE = (1.468.^(1:8)-1)./0.468
diff(log(sumaE))



%% figEpdfhist 
setPlotDefaults 

% summed and local fields comparison
figure
[ah, pos] = tight_subplot(4,1,[.05 .05],[.1 .05],[.1 .05]);

% CHOOSE: pulse size or energy
chplt = 'A2'   % E1 E2 A1 A2
switch chplt % joint local fields: a ~ -2 
  case 'E1', ts = abs(sE); xylim = [.5 100 1e-7 1]; ftl = [1.45 15.3];
  case 'E2', ts = sE.^2; xylim = [1e-1 2e3 1e-7 1]; ftl = [2.1 235];
  case 'A2', ts = sA.^2; xylim = [1e-1 2e3 1e-7 1]; ftl = [1e-2 1e2];  %ts = sA; 
end

axes(ah(3)) % subplot(211)  
[~,edges] = histcounts(log10(ts), 'Normalization', 'pdf'); 
hd=histogram(ts, 10.^edges,'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'EdgeColor', 'k', 'LineWidth', lw);
hold on, axis(xylim), grid on
switch chplt
  case 'E1', ylabel('Probability density ($f_{\mathcal{E}_l}$)', 'Interpreter', 'latex')
  case 'E2', ylabel('Probability density ($f_{\mathcal{E}^2_l}$)', 'Interpreter', 'latex')
  case 'A2', ylabel('Probability density ($f_{A^2_l}$)', 'Interpreter', 'latex')
end
for i=1:nl
    switch chplt
      case 'E1', histogram(sE(:,i), 10.^edges,'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'EdgeColor', 'k', 'EdgeAlpha', .4, 'LineWidth', lw/2)
      case 'E2', histogram(sE(:,i).^2, 10.^edges,'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'EdgeColor', 'k', 'EdgeAlpha', .4, 'LineWidth', lw/2)
      case 'A2', histogram(sA(:,i).^2, 10.^edges,'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'EdgeColor', 'k', 'EdgeAlpha', .4, 'LineWidth', lw/2)
end,end
set(gca, 'xscale', 'log', 'yscale', 'log', 'FontSize', fsz, 'LineWidth', alw)
% fitted line (mnaE2 = 2.22, mnaE1 = 1.45
rgix = 10.^edges>=ftl(1) & 10.^edges<=ftl(2); hc = hd.Values; 
%x = edges(rgix); y = log10(hc(rgix)./diff(10.^[edges(find(rgix,1)-1) edges(rgix)])); 
x = edges(rgix); y = log10(hc(rgix)); pf = polyfit(x,y,1)
hold on, plot(10.^x, 10.^pf(2) * (10.^x).^pf(1), 'k--')
[pks,locs,wdth,prom] = findpeaks(hc,10.^edges(1:end-1)); 10.^diff(log10(locs));

switch chplt % summed field: a ~ -1
    case 'E1', ts = sum(abs(sE),2); ftl = [1.45 44.57]; % same as: ts = abs(sum(sE,2)); 
    case 'E2', ts = sum(sE.^2,2); ftl = [2.1 434]; % mnaE^2
    case 'A2', ts = sum(sA.^2,2); ftl = [1e-2 1e2]; % mnaA^2  % ts = sum(sA,2);
end

axes(ah(4)) % subplot(212) 
[~,edges] = histcounts(log10(ts), 'Normalization', 'pdf');
hd=histogram(ts, 10.^edges, 'DisplayStyle', 'stairs', 'Normalization', 'pdf', 'EdgeColor', 'k', 'LineWidth', lw);
axis(xylim), grid on
switch chplt
  case 'E1', xlabel('Pulse size'), ylabel('Probability density ($f_{\Sigma\mathcal{E}}$)', 'Interpreter', 'latex')
  case 'E2', xlabel('Pulse energy'), ylabel('Probability density ($f_{\Sigma\mathcal{E}^2}$)', 'Interpreter', 'latex')
  case 'A2', xlabel('Subthreshold activity energy'), ylabel('Probability density ($f_{\Sigma A^2}$)', 'Interpreter', 'latex')
end
set(gca, 'xscale', 'log', 'yscale', 'log', 'FontSize', fsz, 'LineWidth', alw) 
% fitted line
rgix = 10.^edges>=ftl(1) & 10.^edges<=ftl(2); hc = hd.Values; 
x = edges(rgix); y = log10(hc(rgix)); pf = polyfit(x,y,1)
%figure, plot(edges(1:end-1),log10(hc)), hold on, plot(x, pf(1)*x+pf(2))
hold on, plot(10.^x, 10.^pf(2) * (10.^x).^pf(1), 'k--')
[pks,locs,wdth,prom] = findpeaks(hc,10.^edges(1:end-1)); 10.^diff(log10(locs));

%%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figAE2pdfhist') % figEpdfhist, figE2pdfhist, figA2pdfhist



%% Visualization: comparing functions of sE at different scales
figure
ts = sum(abs(sE),2); % + sum(abs(L),2);% summed global field, a ~ -1
ts = abs(sE);                          % joint local fields,  a ~ -2
%ts = sum(sE.^2,2);
[~,edges] = histcounts(log10(ts));
subplot(221), histogram(ts), title('E') 
subplot(222), histogram(ts), set(gca, 'yscale', 'log')
subplot(223), histogram(ts, 10.^edges), set(gca, 'xscale', 'log')
subplot(224), histogram(ts, 10.^edges), set(gca, 'xscale', 'log', 'yscale', 'log')

% pulse avalanche size scaling exponent: E      
% ofw ~2.5,2.7;  bfw ~1.5
h = histogram(ts,10.^edges);
%rgix = 10.^edges>=1 & 10.^edges<=75;  
rgix = 10.^edges>=1.5 & 10.^edges<=25; 
% x = (edges(2:end)-0.005) * log(10);
x = edges(rgix); 
y = log10(h.Values(rgix)); %y(isinf(y)) = 1;
pf = polyfit(x,y,1)
figure, plot(edges(1:end-1),log10(h.Values))
hold on, plot(x, pf(1)*x+pf(2))
[pks,locs,wdth,prom] = findpeaks(h.Values,10.^edges(1:end-1))
diff(log10(locs))

%{
figure
x = 1:0.01:80; % x = 1:0.01:1e3;
dfafit = x.^(-alpha) * exp(trcpt); 
hold on, plot(x, dfafit)
ts = sum(abs(Ebf),2); % + sum(abs(L),2);
tsa = abs(ts);
[~,edges] = histcounts(log10(tsa));
subplot(241), histogram(ts), title('Ebf')
subplot(242), histogram(ts), set(gca, 'yscale', 'log')
subplot(243), histogram(tsa, 10.^edges), set(gca, 'xscale', 'log')
subplot(244), histogram(tsa, 10.^edges), set(gca, 'xscale', 'log', 'yscale', 'log')
DL = diff(A,1); B = zeros(size(DL,1),1); for i=1:size(DL,1), j=find(DL(i,:),1,'last'); if isempty(j),j=1; end, B(i)=j; end
Enl = []; for i=1:size(E,1), Enl= [Enl E(i, find(E(i,:),1,'last'))]; end
subplot(2,4,5),  histogram(B), hold on, histogram(L), set(gca, 'yscale', 'log')
subplot(2,4,6), histogram(avlen), set(gca, 'yscale', 'log')
subplot(2,4,7:8), plot(E)
preprint
%}



%% deprecated

% A
sL = A(ni/10+1:end,:); sLbf = Abf(ni/10+1:end,:);
Pa0 = mean(sL==0, 1);
absL = {}; for i=1:nl, absL{i} = abs(sL(sL(:,i)~=0,i)); end
%for i=1:nl, gevd{i} = fitdist(absL{i}, 'GeneralizedExtremeValue'); end
%for i=1:nl, gamd{i} = fitdist(absL{i}, 'Gamma'); end
tabsL = sL(sL~=0);
figure, nr=ceil(sqrt(nl+2));
for i=1:nl+2
    subplot(nr,nr,i)
    if     i == nl+1, histogram(abs(sum(tabsL,2))), title(['abssumL, P(a_{all}=0)=' num2str(mean(Pa0))])
    elseif i == nl+2, histogram(sum(abs(tabsL),2)), title('sumabsL')
    else,  histogram(absL{i}) 
        title(['P(a_' num2str(i) '),  P(a_' num2str(i) '=0)=' num2str(Pa0(i))])
    end, %set(gca, 'yscale', 'log')    
end
preprint

% T
gv = kron(1:nl, ones(1, 0.9*ni))';
sT = T(ni/10+1:end,:); 
rT = reshape(sT, 0.9*ni*nl, 1);
gamd = fitdist(rT,'Gamma', 'By', gv); 
gevd = fitdist(rT,'GeneralizedExtremeValue', 'By', gv); 
figure
for i=1:nl+2
    subplot(nr,nr,i)
    if     i == nl+1, histogram(sum(T,2)), title('sumT')
    elseif i == nl+2, histogram(sum(abs(T),2)), title('sumabsT')
    else,         h = histogram(T(:,i));
        hold on, title(['\theta_' num2str(i) '=' num2str(mdT(i)) '; shape=' num2str(gamd{i}.a) ', scale=' num2str(gamd{i}.b)])
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * gampdf(h.BinEdges, gamd{i}.a, gamd{i}.b))
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gev', h.BinEdges, gevd{i}.k, gevd{i}.sigma, gevd{i}.mu), '--')
    end, %set(gca, 'yscale', 'log')        
end
preprint

% E
sE = E(ni/10+1:end,:); 
Pf = mean(sE~=0, 1); 
absE = {}; for i=1:nl, absE{i} = abs(sE(sE(:,i)~=0,i)); end
for i=1:nl, gevd{i} = fitdist(absE{i}, 'GeneralizedExtremeValue'); end
for i=1:nl, gamd{i} = fitdist(absE{i}, 'Gamma'); end
tabsE = sE(sE~=0);
figure
for i=1:nl+2
    subplot(nr,nr,i)
    if     i == nl+1, histogram(abs(sum(tabsE,2))), title('abssumE')
    elseif i == nl+2, histogram(sum(abs(tabsE),2)), title('sumabsE')
    else,         h = histogram(absE{i});
        hold on, title(['P(e_' num2str(i) '),  Pf=' num2str(Pf(i))])
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gev', h.BinEdges, gevd{i}.k, gevd{i}.sigma, gevd{i}.mu), '--')
        plot(h.BinEdges, sum(h.Values)*h.BinWidth * pdf('gam', h.BinEdges, gamd{i}.a, gamd{i}.b), '-')
    end, %set(gca, 'yscale', 'log')    
end
preprint

