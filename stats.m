ACT=A; A=L; L=V;
%%
ns = floor(ni*.9);
bip = ni*.1;   % burn-in periodACT

sA = A(bip+1:end,:); %sAbf = Abf(bip+1:end,:);
sT = T(bip+1:end,:); 
sE = E(bip+1:end,:); 
gv = kron(1:nl, ones(1, ns))';
rT = reshape(sT, ns*nl, 1);
rE = reshape(sE, ns*nl, 1);

%% Calculate statistics from numerical simulations
% stationary probability of firing 
Pf = mean(sE~=0, 1)

% stationary probability of firing given input 
% P(e~=0) > 0.5 for l>2 because stationarity has not been reached
Pfgi = [mean(sE(:,1)~=0, 1) mean(sE(:,1:end-1)~=0 & sE(:,2:end)~=0, 1) ./ mean(sE(:,1:end-1)~=0, 1)]
Pf ./ [1 Pf(1:end-1)];  % can also be computed as the division by the previous Pf
cumprod(Pfgi) % cumprof(Pfgi)=Pf cumulative product of Pfgi is prob of firing

% mean threshold
mnT = mean(sT, 1), median(sT, 1); th1 = mnT(1); sdT = std(sT, 1); 
ratioTT = mean(sT(:,2:end)) ./ mean(sT(:,1:end-1)) % scaling rate ~3/2

% mean error signals
zlE = {}; absE = {};
for i=1:nl, zlE{i} = sE(sE(:,i)~=0,i); absE{i}=abs(zlE{i}); end
sd0E = std(sE,1);
% mnaE.^2 + sdaE.^2 = sdpE.^2
sdE = cellfun(@std, zlE), sdaE = cellfun(@std, absE)
mnaE = cellfun(@mean, absE);  % meanT as mean of absE in l-1
mdaE = cellfun(@median, absE); 
mnT ./ mnaE


% Activities
Pa0 = mean(sA==0, 1)  % approx stationary probability of not being stimulated after having fired
% sdaA.^2 + mnaA.^2 = sdA.^2
sdA = std(sA, 1) 
mnaA = mean(abs(sA), 1), sdaA = std(abs(sA), 1)  % total mean


% A tilde recomputation: ONLY VALID FOR H=1/2
for i=1:nl
    sAc{i} = sA(:,i); asAc{i} = abs(sAc{i});
    zlA{i} = sA(sA(:,i)~=0, i);   % A for non-zero trials
    azlA{i} = abs(zlA{i}); 
    
    if i==1
        ATwz{1} = sA(randi(ns,[ns,1]), i) + normrnd(0,1,[ns,1]);
        AT{1} = ATwz{1};
        ATr{1} = ATwz{1};        
    else
        nsi = length(absE{i-1});
        nss = ns; % nsi ns
        AT{i} = sA(randi(nsi,[nss,1]), i) + zlE{i-1}(randi(nsi,[nss,1]));
        ATr{i} = sA(randi(nsi,[nss,1]), i) + absE{i-1}(randi(nsi,[nss,1]));
        ATwz{i} = sA(randi(ns,[ns,1]), i) + sE(randi(ns,[ns,1]), i-1);
    end
    aATwz{i} = abs(ATwz{i});
    aAT{i} = abs(AT{i});
    %for j = 1:ns
    %    if i==1,              AT{1}(j) = normrnd(0,1) + sA(j,1);
    %    elseif sE(j, i-1)~=0, AT{i}(j) = sE(j,i-1) + sA(j,i); end
    %end
end

% for non-zero trials
% mnaAnz.^2 + sdaAnz.^2 = sdAnz.^2
sdAnz = cellfun(@std, zlA)  
mnaAnz = cellfun(@mean, azlA), sdaAnz = cellfun(@std, azlA)   
apsA = mnaA ./ (1-Pa0)   % apsA = mnaAnz
atsA = apsA ./ mnT   % for non-zero trials as a fraction of threshold

% A below g^.5*mnaE
prevmnE = num2cell([1 mnaE(1:end-1)]);
mnaAbme = cellfun(@(x,y) mean(x(x <= y)), asAc, prevmnE);
sdaAbme = cellfun(@(x,y) std(x(x <= y)), asAc, prevmnE);
sdAbme = cellfun(@(x,y) std(x(x <= y)), sAc, prevmnE);

% standard deviation of absolute Atilde
sdATwz = cellfun(@std, ATwz)
mnaATwz = cellfun(@mean, aATwz), sdaATwz = cellfun(@std, aATwz) 
% for non-zero E input
% mnaATnz.^2 + sdaATnz.^2 = sdATnz.^2
sdAT = cellfun(@std, AT)
mnaAT = cellfun(@mean, aAT), sdaAT = cellfun(@std, aAT) 

% ATr negative, ATneg
sdATr = cellfun(@std, ATr); % sum(ATr^2) = sum(aAT^2)
mnATr = cellfun(@mean, ATr); 
pATr_neg = cellfun(@(x) mean(x<0), ATr)
st0D = cellfun(@(x) mean(x(x<0).^2).^.5, ATr);
mnD = cellfun(@(x) -mean(x(x<0)), ATr); stD = cellfun(@(x) std(x(x<0)), ATr);


% activity given that the next trial fired 
for i=1:nl
  gnfA{i} = abs(sA([sE(2:end, i)~=0; false], i)); 
  switch i
    case 1,    gnrA{i} = abs(sA([sE(2:end, i)==0; false], i)); 
    otherwise, gnrA{i} = abs(sA([sE(2:end, i)==0; false] & [sE(:, i-1)~=0], i)); 
  end
end
mnaAgnf = cellfun(@mean, gnfA);
mnaAgnr = cellfun(@mean, gnrA);


% Summary
h = 0.5;
[Pf.^-h; sd0E; sdE; mnaE; sdaE; mnT; sdAnz; mnaAnz; sdaAnz; sdAT; mnaAT; sdaAT; sdA; mnaA; sdaA; mnaAgnf]
[Pf.^-h; sd0E; sdE; mnaE; sdaE; mnT; sdAT; mnaAT; sdaAT].*Pf.^h
[sdAnz; mnaAnz; sdaAnz; sdA; mnaA; sdaA; mnaAgnf; sdATr; mnATr; st0D; mnD; stD].*Pf.^h
figure
plot(([Pf.^-0.5; sd0E; sdE; mnaE; sdaE; mnT; sdAnz; mnaAnz; sdaAnz; sdAT; mnaAT; sdaAT; sdA; mnaA; sdaA; mnaAgnf].*Pf.^h)'), grid on, ylim([0 1])
line(repmat([1.5;6.5],1,5),repmat([1/2 1/3 2/3 1/4 3/4], 2, 1), 'LineStyle', ':')


% Partitioning Atilda wrt T_{l-1}
for i=1:nl  
  if i==1, lprA(i) = 0;   
  else
    AbTi = abs(zlA{i}) < mnT(i-1); % pA trials where A_l < T_l-1
    lprA(i) = mean(AbTi);   % fraction of A_l below T_l-1
    mnaAlT(i) = mean(zlA{i}(zlA{i} < mnT(i-1)), 1);
    mnaAuT(i) = mean(zlA{i}(zlA{i} > mnT(i-1)), 1);
  end    
end
g = Pa0;       % actual gain
gb = ([mnT(2:end)./mnT(1:end-1) NaN]).^-2;
ga = 2*(1-normcdf(mnT,0,[1 Pf(1:end-1)].^-.5)); % discrete Gaussian input gain
% lprAT + uprAT + Pa0 = 1
lprAT = (1-Pa0).*lprA; % fraction of Atilda below T_l-1 
uprAT = (1-Pa0).*(1-lprA); % fraction of Atilda above T_l-1


% mnT as the weighted sum of mnaA and mnaE 
(mnaE.*Pa0 + mnaAnz.*(1-Pa0)) ./ mnT  % ~1  
%also (mnaE.*Pa0 + mnaA) = mnT  because AT and A coincide in (0,T)
(mnaE + mnaAnz)/2 ./ mnT   % ~1/0.97
(mnaE + mnaA)/2 ./ mnT   % ~8/9

% Expected values SEM and ratia TT, TE
ratioTT = mnT./[1 mnT(1:end-1)];
diffTE = mnT(2:end)-mnaE(1:end-1);
%ratioTE = [mnT(1)/sqrt(2/pi) mnT(2:end)./mnaE(1:end-1)]; % rTT = rTE(2:end) .* rET(1:end-1)
ratioTE = mnT./[sqrt(2/pi) mnaE(1:end-1)]; % rTT = rTE .* rET % kte
ratioET = [mnaE ./ mnT]; % ket
ratioAT = [mnaA ./ mnT]; 
[ratioTT; ratioTE; ratioET; Pa0; ratioTT.*Pa0.^.5; mnT.*Pf.^.5; ratioAT]'


%% Uncertainty estimates
% Bias corrected and accelerated percentile bootstrap CI estimates (Efron, 1987)
%bci = bootci(100, {@mean, absE{1}}, 'alpha', .05, 'type', 'bca') % for absE1: 1.4467, 1.4479, 1.4499

% autocorrelation time to compute effective sample size
Tact = 2 * 1/w; % 1/w % 2/w-1=1+2*sum1toInf(rho) Straatsma86, Kass98;  
Aact = 2 * 3/4*mnT.^2;
Eact = 1;

%figure, errorbar(1:(nl-1), ratioTE, semE(1:end-1), semT(2:end)), grid on
%figure,plot(ratioTE, skewness(sT(:,2:end)), '*')
for i = 1:nl
    semaA(i) = std(abs(sA(:,i)))./sqrt(length(sA(:,i))./Aact(i));
    sesdA(i) = std(sA(:,i))./sqrt(length(sA(:,i))./Aact(i));
    semaAnz(i) = std(abs(sA(sA(:,i)~=0,i)))./sqrt(length(sA(sA(:,i)~=0,i))./Aact(i));
    semaT(i) = std(sT(:,i))/sqrt(length(sT(:,i))/Tact);
    semaE(i) = std(absE{i})/sqrt(length(absE{i})/Eact);
    if i==1, semrTT(i) = std(sT(:,i))/sqrt(ns);
    else,    semrTT(i) = std(sT(2:end,i)./sT(1:end-1,i-1))/sqrt(ns-1); end
    if i==1, semrTE(i) = std(sT(:,i)/sqrt(2/pi))/sqrt(ns);
    else, isnei0=sE(:,i-1)~=0; sTnei0 = sT(isnei0,i); absEnei0 = abs(sE(isnei0,i-1));
             semrTE(i) = std(sTnei0(2:end) ./ absEnei0(1:end-1)) / sqrt(sum(isnei0)); end
end, semaA, sesdA, semaAnz, semaT, semaE, semrTT, semrTE

% Normal approximation to binomial proportion
nsE = ns*[1 Pf(1:end-1)];
semGain = sqrt(Pfgi.*(1-Pfgi) ./ nsE)  % gain
semFR = sqrt(Pf.*(1-Pf) ./ ns)  % firing rate
semPa0 = sqrt(Pa0.*(1-Pa0) ./ ns) % P(a=0)

% exact binomial CI (computationally intensive)
%binoinv([.05 .95], ns, Pfgi(1))/ns  % Clopper-Pearson?
% binomial fit, MLE?
%[PHAT, PCI] = binofit(sum(sE~=0,1), ns*ones(1,7), .05); [PCI(:,1) PHAT' PCI(:,2)]

%% probability of sign(E+A) ~= sign(E)
sgE=zeros(nl,ni*.9); sgAE=zeros(nl,ni*.9);
for i = 2:nl
    for j = 2:size(sE,1)
        if sE(j-1,i-1) ~= 0
            sgE(i,j) = sign(sE(j-1,i-1));
            sgAE(i,j) = sign(sA(j,i) + sE(j-1,i-1));
        end
    end
    sum(sgE(i,:)~=0)
    mean(sgE(i,sgE(i,:)~=0) ~= sgAE(i,sgAE(i,:)~=0))    
end

%% Utils

% Threshold pdfs
pdlist = {'Nakagami', 'Rician', 'Poisson', 'Lognormal', 'Burr', 'Logistic', 'Weibull'};
pdbic = @(pd) 2*pd.negloglik + pd.NumParameters*log(length(pd.InputData.data));

% using distFit app, Gamma seems to fit best, followed by GEV and Nakagami
pd = fitdist(rT,'Gamma', 'By', gv)   % shape, scale
pd{1}.a * pd{1}.b, pd{1}.negloglik
for i=1:nl, pd{i} = fitdist(absE{i}, 'Gamma'); pd{i}, pd{i}.negloglik, end

% GEV
Tpd = fitdist(rT,'GeneralizedExtremeValue', 'By', gv)
for i=1:nl
    Epd{i} = fitdist(pE{i}, 'GeneralizedExtremeValue'); 
    pdbic(Epd{i})
end
% in theory the specific GEV should be a sth like a Gumbel distribution
gumbel = @(x, sigma, mu) gevpdf(x, 0, sigma, mu); % fix k=0
for i=1:nl, mle(pE{i}, 'pdf', gumbel, 'start', [1, 1]), end

for pdname = pdlist
    pd = fitdist(rT, pdname, 'By', gv); 
    pdbic(pd{1}) 
end


% Estimation of neurofrequency bands
% if we assume an input rate ~ 200Hz (5ms spike + refractory period)
sif = 200; mnL = []; 
for i = 1:nl, mnL = [mnL mean(L==i)]; end
nfbV = mnL * sif % highest layers (6,7) have 4~6Hz, like theta rhythm 4~7
nfbE = Pf * sif


%% Recurrence plot
n = 2e2; % ni
repl = zeros(n,n);

v = cae(1e4:end);
for i=1:n,for j=1:n
  if abs(v(i)-v(j)) < .1
    if v(i)==0 || v(j)==0, repl(i,j)=0.3; %.5
    else,                  repl(i,j)=1; 
    end
  end
end, end
figure, imshow(repl)

%% Poincare map / First recurrence map
n = 1e5; % ni
poimae = zeros(n-1,2); poime = zeros(n-1,2);
for i=1:n-1     
    poimae(i,:) = abs(ca(i:i+1)); %.^2
    poime(i,:) = abs(ce(i:i+1));
end

figure
subplot(231),scatter(poimae(:,1),poimae(:,2)), xlabel('a_t'),ylabel('a_{t+1}'), grid on, box on
%poimh=poimae; poimh(any(poimae==0,2),:)=[];
subplot(232),histogram2(poimae(:,1),poimae(:,2),'DisplayStyle','tile'), xlabel('a_t'),ylabel('a_{t+1}')
subplot(233),histogram2(poimae(:,1),poimae(:,2),'DisplayStyle','bar3'), xlabel('a_t'),ylabel('a_{t+1}')
subplot(234),scatter(poime(:,1),poimae(:,2)), xlabel('e_t'),ylabel('a_{t+1}'), grid on, box on
subplot(235),histogram2(poime(:,1),poimae(:,2),'DisplayStyle','tile'), xlabel('e_t'),ylabel('a_{t+1}')
subplot(236),histogram2(poime(:,1),poimae(:,2),'DisplayStyle','bar3'), xlabel('e_t'),ylabel('a_{t+1}')


%% Transition matrices
% don't forget consecutive values are rescaled by 1/k

jt = eps*sign(randn(length(caij),1)); % jittering eschews artifact at [0,0]

poima1=caij(:,1)+jt; poima2=caij(:,2)+jt;poime1=ceij(:,1)+jt; poime2=ceij(:,2)+jt;
%poima1=abs(poima1); poima2=abs(poima2);poime1=abs(poime1); poime2=abs(poime2);
sga1=sign(poima1); pmaa1 = poima1.*sga1; pmaa2 = poima2.*sga1;
sgr=sign(randn(size(poima1))); pmaa1=pmaa1.*sgr; pmaa2 = pmaa2.*sgr;
sge1=sign(poime1); pmee1 = poime1.*sge1; pmee2 = poime2.*sge1;
sgr=sign(randn(size(poime1))); pmee1=pmee1.*sgr; pmee2 = pmee2.*sgr;
                   pmea1 = pmee1;        pmea2 = poima2.*sge1;
sgr=sign(randn(size(poime1))); pmea1=pmea1.*sgr; pmea2 = pmea2.*sgr;
                   pmae1 = pmaa1;        pmae2 = poime2.*sga1;
sgr=sign(randn(size(poima1))); pmae1=pmae1.*sgr; pmae2 = pmae2.*sgr;

bw=0.05; sz=4; pcsj=[0 .015]; pcsc=[0 .5]; pcsj = [0 .1];
axi = {[-1-bw/2 1+bw/2 -1-bw/2 1+bw/2], [-3 3 -3 3], [-3 3 -1-bw/2 1+bw/2], [-1-bw/2 1+bw/2 -3 3]};
figure, ha = tight_subplot(5,4);
axes(ha(1)), scatter(pmaa1,pmaa2,sz), xlabel('a_t'),ylabel('a_{t+1}/k'), grid on, box on, axis(axi{1})
axes(ha(5)), histogram2(pmaa1,pmaa2,'DisplayStyle','tile','BinWidth',bw,'Normalization','probability'), xlabel('a_t'),ylabel('a_{t+1}/k'), colorbar, caxis(pcsj)
axes(ha(9)), haa=histogram2(pmaa1,pmaa2,'DisplayStyle','bar3','BinWidth',bw,'Normalization','probability','XBinLimits',axi{1}(1:2),'YBinLimits',axi{1}(3:4)); xlabel('a_t'),ylabel('a_{t+1}/k')
axes(ha(2)), scatter(pmee1,pmee2,sz), xlabel('e_tk'),ylabel('e_{t+1}'), grid on, box on, axis(axi{2})
axes(ha(6)), histogram2(pmee1,pmee2,'DisplayStyle','tile','BinWidth',bw,'Normalization','probability'), xlabel('e_tk'),ylabel('e_{t+1}'), colorbar, caxis(pcsj)
axes(ha(10)), hee=histogram2(pmee1,pmee2,'DisplayStyle','bar3','BinWidth',bw,'Normalization','probability'); xlabel('e_tk'),ylabel('e_{t+1}')
axes(ha(3)), scatter(pmea1,pmea2,sz), xlabel('e_tk'),ylabel('a_{t+1}'), grid on, box on, axis(axi{3})
axes(ha(7)), histogram2(pmea1,pmea2,'DisplayStyle','tile','BinWidth',bw,'Normalization','probability'), xlabel('e_tk'),ylabel('a_{t+1}'), colorbar, caxis(pcsj)
axes(ha(11)),hea=histogram2(pmea1,pmea2,'DisplayStyle','bar3','BinWidth',bw,'Normalization','probability'); xlabel('e_tk'),ylabel('a_{t+1}')
axes(ha(4)), scatter(pmae1,pmae2,sz), xlabel('a_t'),ylabel('e_tk'), grid on, box on, axis(axi{4})
axes(ha(8)), histogram2(pmae1,pmae2,'DisplayStyle','tile','BinWidth',bw,'Normalization','probability'), xlabel('a_t'),ylabel('e_tk'), colorbar, caxis(pcsj)
axes(ha(12)),hae=histogram2(pmae1,pmae2,'DisplayStyle','bar3','BinWidth',bw,'Normalization','probability'); xlabel('a_t'),ylabel('e_tk')

hcell = {'haa', 'hee', 'hea', 'hae'}; 
condxpc = {}; bx={}; by={}; hcv = {}; condabxpc={};
for j=1:4, h=eval(hcell{j});
    % Conditional expectations 
    bx{j} = diff(h.XBinEdges)/2 + h.XBinEdges(1:end-1);
    by{j} = diff(h.YBinEdges)/2 + h.YBinEdges(1:end-1);
    for i=1:h.NumBins(1)        
        condxpc{j}(i) = h.Values(i,:) * by{j}' / sum(h.Values(i,:)); 
        condabxpc{j}(i) = h.Values(i,:) * abs(by{j})' / sum(h.Values(i,:)); 
    end
    % recompute 3D hist with cond pdsif   
    hcv{j} = h.Values;
    for i=1:h.NumBins(1)
        hcv{j}(i,:) = hcv{j}(i,:) / sum(hcv{j}(i,:)); 
    end    
end
% replot
colcon = (exp(8*(1:-.1:0))-1)/(exp(8)-1);
for i=1:4
    axes(ha(i)), hold on, 
    plot(bx{i},condabxpc{i}, 'LineWidth',3), plot(bx{i},condxpc{i}, 'LineWidth',3)
    axes(ha(12+i))
    %bh=bar3(rot90(hcv{i},3)); colormap(flipud(gray))% surfc(hcv{i})
    %for k = 1:length(bh), bh(k).CData = bh(k).ZData; bh(k).FaceColor = 'interp'; end
    [mx,my]=meshgrid(bx{i},by{i}); 
    surf(mx,my,hcv{i}'), %caxis([0 1])
    ca=gca; ca.XLim=ha(i).XLim; ca.YLim=ha(i).YLim;
    axes(ha(16+i))
    colormap(flipud( (exp(8*gray)-1)/(exp(8)-1))) 
    %colormap(repmat(1000.^linspace(0,-1,64)',1,3))
    %imagesc(rot90(hcv{i},3))
    contourf(mx,my,hcv{i}',colcon) 
    colorbar, caxis(pcsc), hold on   
end

%%
%figure, plotmatrix([pmaa1 pmaa2 pmee1 pmee2])
figure
for i=1:4, subplot(2,2,i)
    [mx,my]=meshgrid(bx{i},by{i});
    h=meshc(mx,my,hcv{i}'); h(2).ContourZLevel=1;
    %h=surfc(mx,my,hcv{i}'); %colormap(flipud(gray))
    %h=contour(mx,my,hcv{i}');
    % ksdensity
    ca=gca; ca.ZLim = [0 1];
    colorbar
end

for i=1:4, 
    max(condabxpc{i}), bx{i}(find(max(condabxpc{i})==condabxpc{i}))
    min(condabxpc{i}), bx{i}(find(min(condabxpc{i})==condabxpc{i})) 
end
[srt,idx]=sort(hcv{1}(:,21)), bx{1}(idx)
[srt,idx]=sort(hcv{4}(:,58)), bx{4}(idx)
    
figure, for i=1:4, h=eval(hcell{i}); subplot(2,2,i),plot(h.Values), set(gca,'YScale','log'),end
figure, for i=1:4, subplot(2,2,i),plot(hcv{i}), end

%% Appendix A: modeling A~ and estimating Theta

% a, theta pdf
m=0; s=1;
nd = makedist('normal', 'mu', m, 'sigma', s);
%ud = makedist('uniform', 'lower', 0, 'upper', th1);
%bnd = nd.truncate(-th1, th1);
%utnd = nd.truncate(th1, inf); ltnd = nd.truncate(-inf, -th1);
% truncated gaussian
hnd = makedist('halfnormal', 'mu', m, 'sigma', s);
%bhnd = hnd.truncate(0, th1);
%thnd = hnd.truncate(th1, inf);

% L1: lim t->inf, steady state
% convolve uniform with gaussian: ucg is the limit distribution only for l=L->inf
syms x y pi t x1 x2 
sf = 1/(2*t*sqrt(2*pi)) * int(exp(-(y-x)^2/2), y, [-t, t]); % convolution N*U
sf1 = 2 * sf;   % take as reflected over 0
sf2 = int(sf1, x, [0, t]);
sf3 = int(sf1, x, [t, inf]);
ucg_pdf = matlabFunction(sf1); 
ucg_cdf = matlabFunction(sf2); 
ucg_bmn = matlabFunction(mean(sf2)); 
ucg_tmn = matlabFunction(mean(sf3)); 


% analytical formula for (constant) thresholds in the steady state
iste = 0;
% gaussian input: P(fire) * dist2a = P(rest) * dist2a 
gauss_t1of = @(x) abs((1-hnd.cdf(x)) * (mean(hnd.truncate(x, inf)) - x) - hnd.cdf(x) * (x - mean(hnd.truncate(0, x))));
% truncated gaussian random walk input: P(fire) * dist2a = P(rest) * dist2a + P(rest|prevfire) * x 
syms x y th xb xl
N = symfun(1/sqrt(2*pi) * exp(-x^2/2), x);
nl = 1;
ni = 4;
for l = 1:nl
    switch l
        case 1,    I = N;   % truncated N, for x<theta
        otherwise, I = Pat(l-1);       
    end
    
    fati(l,1) = I;  
    Fati(l,1) = int(fati(l,1), [-th, th]);
    %fati(l,2) = int(subs(fati(l,1), x-y) * N(y) / Fati(l,1), y, [-th, th]); % truncated N*fapT1, for x<theta
    %Fati(l,2) = int(fati(l,2), [-th, th]);
    
    for i = 2:ni
        if iste && i>2
            te = taylor(subs(fati(l,i-1), x-y) * N(y) / Fati(l,i-1));
            fati(l,i) = int(te, y, [-th, th]); 
        else
            % with convolution
            fatin = fati(l,i-1) / Fati(l,i-1);
            fati(l,i) = int(N(y) * subs(fatin, x-y), y, [-th, th]); % truncated N*fapT1, for x<theta
            % with Fourier transform (characteristic function)
            %fati(l,i) = ifourier(fourier(N(x)) * fourier(fati(l,i-1) * rectangularPulse(-th, th, x) / Fati(l,i-1)));
        end
        if iste && i>1
            te = taylor(fati(l,2));
            Fati(l,i) = int(te, [-th, th]);
        else
            Fati(l,i) = int(fati(l,i), [-th, th]);
        end        
    end
    
    % P(a_tilda) up to ti=2
    %Pat(l) = sum(Fati(l,:) .* fati(l,:)) / ni;     
    Pat(l) = Fati(l,1)*fati(l,1) + (1-Fati(l,1)) * (Fati(l,2)*fati(l,2) + 0);
    Pat(l) = Pat(l) / (Fati(l,1) + (1-Fati(l,1))*Fati(l,2)); % normalize
    %ezsurf(Pat)
    % P(a_tilda) up to ti=2
    Pat(l) = 0; Z = 0;
    for j = 2:-1:1
        Pat(l) = Pat(l)*(1-Fati(l,j)) + Fati(l,j)*fati(l,j); 
        Z = Z*(1-Fati(l,j)) + Fati(l,j); 
    end, Pat(l) = Pat(l) / Z; % normalize

    % average probabilities of firing and resting over all trials, up to ti=3
    Pf = 0; Pr = 0;
    for j = 4:-1:1
       Pf = Pf*(1-Fati(l,j)) + Fati(l,j)*(1-Fati(l,j)); 
       Pr = Pr*(1-Fati(l,j)) + Fati(l,j)^2; 
    end
    
    double(subs(Pf,th,mnT(1)))
    double(subs(Pr,th,mnT(1))) 
end
   
% probability of resting (not firing) at each iteration
vpa(subs(Fati, th, 1))
    
% survival function
Sati = 1 - Fati;

% estimate theta through P(a_tilda)
% compute P(a_tilda) body and tails (across theta)
%tail = int(Pat(l), x, [th,inf]) * (int(Pat(l)*x, x, [th, inf]) / int(Pat(l), x, [th,inf]) - th);
%body = int(Pat(l), x, [0,th]) * (th - int(Pat(l)*x, x, [0, th]) / int(Pat(l), x, [0, th]));
tail = int(Pat(l) * (x - th), x, [th, inf]);  %100]);
body = int(Pat(l) * (th-x), x, [0, th]);
rest_eqn = matlabFunction(tail - body);
tgrw_t1of = @(thn) abs(rest_eqn(thn));

% estimate theta through P(a)*I
%tail; body; % <-approximation. P given that there was no pulse in the previous trial: should be slightly higher. arduous to compute
%rest_eqn = matlabFunction((1-Pf)*(tail - body));
%Pfn = matlabFunction(Pf);
%fire_eqn = @(thn) Pfn(thn) * ((1-hnd.cdf(thn)) * (mean(hnd.truncate(thn, inf)) - thn) - hnd.cdf(thn) * (thn - mean(hnd.truncate(0, thn))));
%tgrw_t1of = @(thn) abs(rest_eqn(thn) + fire_eqn(thn));

% numeric optimization of theta
fmbo = optimset('Display', 'iter', 'TolX', 1e-6, 'MaxFunEvals', 500, 'MaxIter', 500);
fminbnd(tgrw_t1of, eps, 10, fmbo)


%% Appendix B.1: Pfiring and B.2: Hitting times

% Hitting time estimated from simulations
dl0 = diff(sA(:,1)==0)';
if dl0(find(dl0~=0, 1, 'last'))==1,      dl0 = [1 dl0 -1]; 
elseif dl0(find(dl0~=0, 1, 'last'))==-1, dl0 = [1 dl0]; 
end
avlen = find(dl0 == -1) - find(dl0 == 1);
h = histogram(avlen);
h.BinCounts * (1:13)' / sum(h.BinCounts)

% P of firing given input at l=1, analytically, using Wald's identity
% Levy distribution: a=.5, b=1 
%levy1 = makedist('Stable','alpha',0.5, 'beta',1, 'gam', th, 'delta',0); 
levy1 = @(x,g) sqrt(g/(2*pi)) * exp(-g./(2*x)) .* x.^(-3/2);
levy1int = @(x,g) integral(@(y) levy1(y,g), x-1, x); %integral(@(x) levy1(x,th1^2), 0, 10^6)
levy1exc = @(x,g) levy1int(x,g) / integral(@(y) levy1(y,g), x-1, inf);

%  first exit time
Fet = @(x,g) 1 - vpa(subs(Fethn,{th, n, t},{g, 10, x}));
fet = @(x,g) vpa(subs(fethn,{th, n, t},{g, 10, x}));

% 
fB = @(x, t) 1/sqrt(2*pi*t) * exp(-x.^2/(2*t));
Phi10 = @(th) 1 - integral(@(x) fB(x,1), -th, th);                   % lower bound for Phi1
Phi1Bt = @(th, t) 1 - integral(@(x) fB(x,t), -th, th);               % upper bound for Phi1
fBU = @(x, th) integral(@(y) fB(y,1) .* heaviside(+th+x-y).*heaviside(+th-x+y)./(2*th), -th, th);
Phi1U = @(th) integral(@(x) fBU(x, th), -th, th, 'ArrayValued', 1);  % better upper bound for Phi1
% Phi10 for i=1, PhiU1 otherwise                                    % even better upper bound 

% Ornstein-Uhlenbeck process pdf
fOU = @(x, t, b, s) sqrt(b./(pi*s.^2*(1-exp(-2*b*t)))) * exp(-b./s.^2 * x.^2./(1-exp(-2*b*t)));
Phi1OU = @(th, t, b) 1 - integral(@(x) fOU(x,t,b,1), -th, th); 

% GPD
sympref('HeavisideAtOrigin', 0.5)
xi=-1/3;
sgm = @(xi, th, k) -xi*th.*(2*k-1); 
%maxA = @(xi, th, k)  th - sgm(xi, th, k) ./ xi;  % 2*k*th
fGpd = @(x, mu, sigma, xi) heaviside(x-mu) .* (1 + xi*(x - mu)./sigma).^(-(xi+1)/xi) ./ sigma .* heaviside((mu-sigma/xi)-x); 
%SGpd = @(x, mu, sigma, xi) (1 + xi*(x-mu)./sigma).^(-1/xi);
SGpd = @(x, mu, sigma, xi) (1 + xi*(x-mu)./sigma).^(-1/xi) .* heaviside(x-mu) + heaviside(mu-x);
FGpd = @(x, mu, sigma, xi) 1-SGpd(x,mu,sigma,xi);
QGpd = @(p, mu, sigma, xi) mu + sigma*(p.^-xi - 1)/xi; % if xi=0, x = mu - sigma*log(p)  
SGpd(3-2^(2/3),1,2/3,-1/3)
gpdmn = @(mu,sigma,xi) mu + sigma./(1-xi); gpdmn(1,2/3,-1/3)
gpdmd = @(mu,sigma,xi) mu + sigma.*(2.^xi-1)./xi; gpdmd(1,2/3,-1/3)

% GPD noise, activity suppressed every trial
PhiGpd0 = 1 - integral(@(x) fGpd(x,1,2/3,-1/3), 0, 3/2);   % .75^3


%% Appendix B.2: Symbolic solver, first exit times

%assume(x > -th & x < th)
S = solve(int(Pat, x, [0, th]) == int(Pat, x, [th, inf]));
vpa(S)

% comparison with simulations
double(int(subs(Pat,th,0.8373),x,[-0.8373,0.8373]))  % prob of firing at layer 1, 0.415


% Reflected P(a12) distribution
syms x pi t xi y sigma z w
N = symfun(1/sqrt(2*pi) * exp(-x^2/2), x); % normal
GP = (xi*(x-t)/sigma + 1)^(-(xi+1)/xi) ; % generalized Pareto
gp = (xi*z + 1)^(-(xi+1)/xi);

gpc = int(subs(gp,y) * subs(gp,z-y), y, [-inf,inf]);
gpct = taylor(gpc);

1/(2*t*sqrt(2*pi)) * int(exp(-(y-x)^2/2), y, [-t, t]); % convolution N*U


% distribution GPD*Unif
syms x pi t xi y sigma z w
z = (x - t) / sigma;
gp = 1/sigma * (xi*z + 1)^(-(xi+1)/xi) * heaviside(x-t); % same as GP
sigma=2/3; xi= -1/3; t=1;
gp = subs(gp, x, abs(x))  % ezplot(gp)
fgu = int(subs(gp,x,y) * heaviside(th+x-y) * heaviside(th-x+y)/(2*th), y, [-th,th]);
pfgu = 1-int(subs(fgu,th,3/2), x, [-3/2,3/2])  % .75^4


% First exit time from strip for Brownian motion
syms th pi t k n x y H
Fethn10 = 1/2 * 2/sqrt(pi) * symsum(int(exp(-x^2), [th*(4*k-1)/sqrt(2*t), th*(4*k+1)/sqrt(2*t)]) - int(exp(-x^2), [th*(4*k-3)/sqrt(2*t), th*(4*k-1)/sqrt(2*t)]), k, -n, n);
Fethn11 = -1 + 2/sqrt(pi) * symsum(int(exp(-x^2), [th*(4*k-1)/sqrt(2*t), th*(4*k+1)/sqrt(2*t)]), k, -n, n);
fethn10 = th / (2*pi*t^3)^.5 * symsum((4*k+1)*exp(-th^2 * (4*k+1)^2 / (2*t)) - (4*k-1)*exp(-th^2 * (4*k-1)^2 / (2*t)), k, -n, n);
% B1 and B3
fethn1 = th * sqrt(2/(pi*t^3)) * symsum((-1)^k * (2*k+1) * exp(-th^2 * (2*k+1)^2 / (2*t)), k, 0, n); 
Fethn1 = int(subs(fethn1,t,y), y, [0, t]);
Sethn1 = 1 - Fethn1;
fethn2 = pi/(2*th^2) * symsum((-1)^k * (2*k+1) * exp(-pi^2 * (2*k+1)^2 * t / (8*th^2)), k, 0, n); 
Fethn2 = int(subs(fethn2,t,y), y, [0, t]);
Sethn2 = 1 - Fethn2;
% IT WORKS! but the second converges faster
% the mean is correct (th1^2)1
fethn1mn = int(subs(fethn1,t,y)*y, y, [0, inf]); vpa(subs(fethn1mn,{th, n},{0.8464, 10}))
fethn2mn = int(subs(fethn2,t,y)*y, y, [0, inf]); vpa(subs(fethn2mn,{th, n},{0.8464, 10}))
th1 = 0.8464;
vpa(subs(fethn1,{th, n, t},{th1, 10, 1}))
vpa(subs(fethn2,{th, n, t},{th1, 10, 1}))

% ksi to k relationship
figure
ezplot('1/2+sqrt(.25+1/(1-ksi))'), grid on

fourier(gp, z, w)

%% Appendix B.3

% numeric optimization of OU drift parameter
bof = @(b) abs(PhiCalc(@(th, t) Phi1OU(th, t, b), th1, 30) - Pa0(1));
fmbo = optimset('Display', 'iter', 'TolX', 1e-6, 'MaxFunEvals', 500, 'MaxIter', 500);
fminbnd(bof, eps, 10, fmbo)
% drift is 1389721134
% OUP expected hitting time: ouht = @(th,x) 1/2 * sum(2.^((1:x)/2) .* th.^(1:x) .* gamma((1:x)/2) ./ factorial(1:x));


%% Appendix D: Monte Carlo checks

ni = 10^6;
nl = 3;
w = 0.01;
gp1 = makedist('GeneralizedPareto', 'k',-1/3, 'sigma',2/3, 'theta',1); %      x^2

mcA = zeros(ni, nl); mcT = zeros(ni, nl); mcE = zeros(ni, nl); mcAE = zeros(ni, nl);  
mca = zeros(1, nl); wmnt = zeros(1, nl); mce = zeros(1, nl); mcae = zeros(1, nl);
for i = 1:ni
    %mca(1) = mca(1) + normrnd(0, 1);  
    %mca(1) = mca(1) + normrnd(1, .1) * (2*randi(2)-3);  
    mca(1) = mca(1) + gp1.random *(2*randi(2)-3);
    
    mcae(1) = mca(1);
    wmnt(1) = wmnt(1) * (1-w) + w * abs(mca(1));
    if abs(mca(1)) > wmnt(1)
        mca(2) = mca(2) + mca(1); mcae(2) = mca(2);
        wmnt(2) = wmnt(2) * (1-w) + w * abs(mca(2));    
        if abs(mca(2)) > wmnt(2)
            mca(3) = mca(3) + mca(2); mcae(3) = mca(3);
            wmnt(3) = wmnt(3) * (1-w) + w * abs(mca(3)); 
            if abs(mca(3)) > wmnt(3)
                %
                mce(3) = mca(3); 
                mca(3) = 0;                
            end
            mce(2) = mca(2); 
            mca(2) = 0;
        end
        mce(1) = mca(1); 
        mca(1) = 0;
    end    
    mcA(i,:) = mca; mcAE(i,:) = mcae; 
    mcT(i,:) = wmnt; mcE(i,:) = mce;
end

gc=mean(mcA(ni/10:end,:)==0) % l1w.1: prob of firing 0.4294
meanT = mean(mcT(ni/10:end,:))    % l1w.1: mean threshold 0.8489
meanT./[1 meanT(1:2)]
medianT = median(mcT(ni/10:end,:))  % l1w.1: median threshold 0.8396
mnabsA = mean(abs(mcA(ni/10:end,:)))
for i=1:nl, mcAi=mcA(ni/10:end,i); mnabsAnz(i) = mean(abs(mcAi(mcAi~=0))); end, mnabsAnz
mnabsAE = mean(abs(mcAE(ni/10:end,:))) 


figure
for l=1:nl
    smA = mcA(:,l); smA(smA==0) = [];
    smT = mcT(:,l); smT(smT==0) = [];
    smE = mcE(:,l); smE(smE==0) = [];
    smAE = mcAE(:,l); smAE(smAE==0) = [];
    subplot(nl, 4, 1 + 4*(l-1))
    histogram(smAE), xlim([-4*1.5^l 4*1.5^l])
    by = get(gca, 'YLim'); line([meanT(1:nl); meanT(1:nl)], [zeros(1,nl); by(2)*ones(1,nl)], 'LineStyle', '--')
    grid on
    subplot(nl, 4, 2 + 4*(l-1))
    histogram(smA), %xlim([-4*1.5^l 4*1.5^l])
    by = get(gca, 'YLim'); line([meanT(1:nl); meanT(1:nl)], [zeros(1,nl); by(2)*ones(1,nl)], 'LineStyle', '--')
    grid on
    subplot(nl, 4, 3 + 4*(l-1))
    histogram(smT), %xlim([-4*1.5^l 4*1.5^l])
    by = get(gca, 'YLim'); line([meanT(1:nl); meanT(1:nl)], [zeros(1,nl); by(2)*ones(1,nl)], 'LineStyle', '--')
    grid on
    subplot(nl, 4, 4 + 4*(l-1))
    histogram(smE), %xlim([-4*1.5^l 4*1.5^l])
    by = get(gca, 'YLim'); line([meanT(1:nl); meanT(1:nl)], [zeros(1,nl); by(2)*ones(1,nl)], 'LineStyle', '--')
    grid on
end


%% cascade of truncated energy pdfs with GPD driving input

gp1 = makedist('GeneralizedPareto', 'k',-1/3, 'sigma',2/3, 'theta',1); 
mea0 = gp1.mean;
sde0 = sqrt(mea0^2+gp1.std^2);

nl = 4; ni = 1.5e5;
cv = zeros(5,ni); ca = zeros(5,ni); cae = zeros(5,ni); ce = zeros(5,ni);
cv(1,:) = gp1.random(ni,1) .* (2*randi(2,ni,1)-3);
ce(1,:) = cv(1,:);  
for j=1:nl
  % shift and rescale?
  nz = ce(j,:)~=0; cenz = ce(j, nz);
  %shift = mnae(j) - 1    % delta
  resca = std(cenz) / sde0; % ktt ; mean(nz)^-.5 for BM  
  %cer = cenz - shift * sign(cenz); 
  cer = cer / resca; 
  
  for i=1:ni
    cv(j+1,i) = cv(j,i) + gp1.random *(2*randi(2)-3);        
    %
    %cae(j,i) = ca(j,i) + gp1.random *(2*randi(2)-3);
    cae(j,i) = ca(j,randi(ni,1)) + cer(randi(length(cer),1));
    %
    cmt(j) = mean(abs(cae(j,1:i))); % gm
    switch abs(cae(j,i)) <= cmt(j)
        case true,  ca(j,i) = cae(j,i); ce(j+1,i) = 0; 
        case false, ce(j+1,i) = cae(j,i); ca(j,i) = 0;
    end
  end
end
% std(abs(ce')).^2 + mean(abs(ce')).^2 = std(ce').^2
cmt
std(ce')
mean(ce'~=0)

figure 
for i=1:nl
    subplot(nl,4,4*i-3), histogram(cv(i,:))
    grid on, xlim([-10 10]), title(num2str(gp1.k))
    subplot(nl,4,4*i-2), histogram(cae(i,:))
    grid on, xlim([-6 6]), title(mean(abs(cae(i,:))))
    subplot(nl,4,4*i-1), nz=ca(i,:)~=0; histogram(ca(i,nz))
    grid on, xlim([-3 3]), title([num2str(mean(~nz)) ' ' num2str(mean(abs(ca(i,nz))))])
    subplot(nl,4,4*i), nz=ce(i,:)~=0; histogram(ce(i,nz))
    grid on, xlim([-6 6]), title(mean(abs(ce(i,nz))))
end

%% Appendix D: self-feeding cascade with rescaling

%w = 0.01, H=0.5
%
hr = 0.025;
for j=1:8
for H = 0:hr:1, [j H]
%fstr='sim_sfl/v0/rs_H%.2f_ni1e5wi1e5_w%1.0e.mat'; fname=sprintf(fstr,H,w);
fstr='sim_sfl/rs_H%.2f_ni1e5_wi%1.0e_mean.mat'; fname=sprintf(fstr,H,wi);
%fstr='sim_sfl/rs_H%.2f_ni1e5wi1e5_mw%1.0e_w%1.0e.mat'; fname=sprintf(fstr,H,mw,w);
if exist(fname,'file'), load(fname, '-regexp','^(?!H$)\w'), fname
else, load(sprintf(fstr,H-hr,wi),'-regexp','^(?!H$)\w'), H-hr 
end   
%}

nw = 1e5; 
wi = 1e4;  % ni, ni/10, ni/100
mw = 1;    % sgd weighting window length
w  = 1e-2; % equivalent wi = 100 ?

% jumble and resize
ce=ce(randi(ni,1,nw)); ca=ca(randi(ni,1,nw)); cae=cae(randi(ni,1,nw)); ct=ct(randi(ni,1,nw));
wmne=wmne(randi(ni,1,nw)); wmng=wmng(randi(ni,1,nw)); wmnt=wmnt(randi(ni,1,nw)); wmxe=wmxe(randi(ni,1,nw)); wsda=wsda(randi(ni,1,nw));
ni = length(ce);
ceij=zeros(ni,2); caij=zeros(ni,2); 

for i=1:ni 
    % adjust skewness
    if mod(i,ni/100)==1  %ni/100
      if 0 % H==0.5 % 0, remove skewness
        %ca=ca.*(2*randi(2,1,ni)-3);  
        ce=ce.*(2*randi(2,1,ni)-3); %cae=cae.*(2*randi(2,1,nw)-3);
      else  % H is probability of having same sign
        ca=abs(ca);
        ce=abs(ce).*sign(double(rand(1,ni)<H)-.5);
      end
    end

    % jumble sample order with replacement
    %ri = randi(ni,1,wi); % yoke or pick independently
    %poa = ca(ri);   % ri, randi(ni,1,wi)
    %poe = ce(ri);   % ri, randi(ni,1,wi)
    %poae = cae(ri); % ri, randi(ni,1,wi)
    % window shifting
    csi = circshift(1:ni,wi-i+1); 
    poa = ca(csi(1:wi)); poe = ce(csi(1:wi)); poae = cae(csi(1:wi)); 
    csi = circshift(1:ni,mw-i+1); 
    pot = ct(csi(1:mw)); 
    
    % rescaling factor and energy addition
    resca(i) = mean(abs(poe)); % mean(abs(ce)), std(ce), 
    % subthreshold activity + pulse
    ceij(i,1) = poe(randi(wi,1));
    caij(i,1) = poa(randi(wi,1)); % a,e randomized
    %caij(i,1) = poa(1);    % a serially linked  % not needed anymore because of symmetry?
    cae(i) = caij(i,1) + ceij(i,1)/resca(i);  
    % threshold
    if i<2, ct(1) = tht/(gst^.5*mae);
    else,   ct(i) = ct(i-1)*(1-w) + w*abs(cae(i-1)); end 
    
    % threshold update rule 
    % fixT:std(poe)*tht, free:wmnt(i)  
  wmnt(i) = mean(abs(poae)); % meanabs, medianabs, std, 1   
  %wmnt(i) = mean(abs(pot)); % stogradesc, g closer to g* using poae than using pot?
    %wmnt(i) = tht*std(ce);    % peg T to sdE
    if abs(cae(i)) <= wmnt(i)  % rest
        ca(i)=cae(i);         
        ceij(i,2) = 0; % ce(i)=0; %not modeling this
    else                       % fire    
        ce(i)=cae(i); ca(i)=0;         
        ceij(i,2) = ce(i);
    end
    
    % bookkeeping
    wsda(i) = sqrt(mean(poa.^2)); % also std(poa) if mean(poa)=0
    wmng(i) = mean(poa==0);
    wmne(i) = mean(abs(poe));
    wmxe(i) = max(abs(poe));
    caij(i,2) = ca(i);
end
save cer4.mat ca cae ct ce wmnt wmxe wmng wi wmne wsda w H ni mw ceij caij

%
movefile('cer4.mat',fname,'f')
end
end
%}


%%
cg = mean(wmng), cg2=mean(ca==0)
mnace = mean(abs(ce)) % rescaled to 1
sdace = std(abs(ce))
sdce = sqrt(mean(ce.^2))  % sqrt(mnace^2+sdace^2); std(ce)
mnaca = mean(abs(ca)), sdaca = std(abs(ca)), sdca = mean(ca.^2)^.5 % sdca=std(ca)
mean(abs(ca(ca~=0)))/sdce
mnacae = mean(abs(cae)), mdcae = median(abs(cae))
cket = mnace/mnacae 
ckte = mnacae, ck = mnace            % th mean, resca mean
%ckte = mnacae*sdce/mnace, ck = sdce  % th mean, resca std
%ckte=mdacae, cket=mnace/mdacae, ck=mnace  % th median, resca mean

mnacae / sdce % why 2/3 ?
sdca / sdce % std(ca); sqrt(mean(ca.^2))
%[rho,pval] = corr([wmng;wmne;wmnt]')
%[rho,pval] = partialcorr([wmng;wmne;wmnt]')

de = 0.46557123187676802665673122521994;
d = [gst, mae*sdce, tht*sdce]; % [de, de+1, 1]; %d = [0.5, 1, 1] 
e = [de, de^-.5, tht*sdce]; %[0.5, sqrt(2), 1]; 

figure
subplot(321), histogram(cae)
grid on, xlim([-4 4]), title(mnacae)
subplot(322), plot(cae(ni*.999:end)), grid on
subplot(323), histogram(ca(ca~=0))
grid on, xlim([-4 4]), title([num2str(cg2) ' '])
subplot(324), plot(ca(ni*.999:end)), grid on
subplot(325), h = histogram(ce);
grid on, xlim([-4 4]), title(mnace), hold on, gpf = fitdist(abs(ce)','gp', 'theta', min(abs(ce))-eps);
plot(h.BinEdges, 1/2*sum(h.Values)*h.BinWidth * pdf('gp', h.BinEdges, gpf.k, gpf.sigma, gpf.theta), '--')
subplot(326), plot(ce(ni*.999:end)), grid on
%
pma=max([wmng;wmne;wmnt]'); pmi=min([wmng;wmne;wmnt]');
axl = [min([d(1),e(1),pmi(1)]) max([d(1),e(1),pma(1)]) min([d(2),pmi(2)]) max([d(2),pma(2)]) min([d(3),pmi(3)]) max([d(3),pma(3)])];
radius = mean(axl(2:2:6)-axl(1:2:6))/2;
axl = axl + ones(1,6)*radius.*[-1 1 -1 1 -1 1];
pma=max([wsda;wmnt;wmne]'/sdce); pmi=min([wsda;wmnt;wmne]'/sdce);
axn = [min([sda,pmi(1)]) max([sda,pma(1)]) min([tht,pmi(2)]) max([tht,pma(2)]) min([mae,pmi(3)]) max([mae,pma(3)])];
radius = mean(axn(2:2:6)-axn(1:2:6))/2;
axn = axn + ones(1,6)*radius.*[-1 1 -1 1 -1 1];
%
figure
subplot(421), plot(wmng), grid on, line([1 ni],[gst gst])
subplot(422), plot(wmng.*(wmng+(wsda/sdce).^2)+(wmnt/sdce).^4), grid on, line([1 ni],[gst gst])
subplot(423), plot(wsda/sdce), grid on, line([1 ni],[sda sda])
subplot(424), plot(wmnt/sdce.*wmng.^.5.*(1-wmng).^.5), grid on, line([1 ni],[sda sda])
subplot(425), plot(wmnt/sdce), grid on, line([1 ni],[tht tht])
subplot(426), plot([(wmng.*wmne/sdce).^.5;wsda./(wmng.*(1-wmng)).^.5/sdce]'), grid on, line([1 ni],[tht tht])
subplot(427), plot(wmne/sdce), grid on, line([1 ni],[mae mae])
subplot(428), plot([((wmnt/sdce).^2+wmng.^2-(wsda/sdce).^2.*(1-wmng))./(2*wmng.*wmnt/sdce); wmng.^-.5/sdce]'), grid on, line([1 ni],[mae mae])
%figure
%subplot(211), plot(wsda./wmnt), grid on
%subplot(212), plot(wmxe/sdce), grid on

%figure
%plotmatrix([wmng;wsda/sdce;wmnt/sdce;wmne/sdce]')
%figure, plotmatrix([ca/sdce;cae/sdce;ce/sdce]')
figure
subplot(4,4,1), histogram(wmng)
subplot(4,4,2), histogram2(wsda/sdce,wmng,50,'DisplayStyle','tile')
subplot(4,4,3), histogram2(wmnt/sdce,wmng,50,'DisplayStyle','tile')
subplot(4,4,4), histogram2(wmne/sdce,wmng,50,'DisplayStyle','tile')
subplot(4,4,6), histogram(wsda/sdce)
subplot(4,4,7), histogram2(wmnt/sdce,wsda/sdce,50,'DisplayStyle','tile')
subplot(4,4,8), histogram2(wmne/sdce,wsda/sdce,50,'DisplayStyle','tile')
subplot(4,4,11), histogram(wmnt/sdce)
subplot(4,4,12), histogram2(wmne/sdce,wmnt/sdce,50,'DisplayStyle','tile')
subplot(4,4,16), histogram(wmne/sdce)
subplot(4,4,[9 10 13 14])
plot3(wsda/sdce,wmnt/sdce,wmne/sdce), grid on, hold on, axis(axn)
plot3(axn(2)*ones(1,ni),wmnt/sdce,wmne/sdce)
plot3(wsda/sdce,axn(4)*ones(1,ni),wmne/sdce)
plot3(wsda/sdce,wmnt/sdce,axn(5)*ones(1,ni))
plot3(sda,tht,mae,'ko','markersize',radius*1e3), plot3(sda,tht,mae,'k*')
plot3(axn(2),tht,mae,'k*'), plot3(sda,axn(4),mae,'k*'),plot3(sda,tht,axn(5),'k*')
%%
figure
plot3(wmng,wmne,wmnt), grid on, hold on, axis(axl)
plot3(axl(2)*ones(1,ni),wmne,wmnt)
plot3(wmng,axl(4)*ones(1,ni),wmnt)
plot3(wmng,wmne,axl(5)*ones(1,ni))
plot3(d(1),d(2),d(3),'ko','markersize',radius*1e3), plot3(d(1),d(2),d(3),'k*')
plot3(axl(2),d(2),d(3),'k*'), plot3(d(1),axl(4),d(3),'k*'),plot3(d(1),d(2),axl(5),'k*')
plot3(e(1),e(2),e(3),'go','markersize',radius*1e3), plot3(e(1),e(2),e(3),'g*')
plot3(axl(2),e(2),e(3),'g*'), plot3(e(1),axl(4),e(3),'g*'),plot3(e(1),e(2),axl(5),'g*')
figure
tt=floor(ni/9); for i=1:9
  subplot(3,3,i), i0 = tt*(i-1)+1; 
  plot3(wmng(i0:i0+tt),wmne(i0:i0+tt),wmnt(i0:i0+tt))
  grid on, hold on, axis(axl)
  plot3(d(1),d(2),d(3),'ko', 'markersize', 10),plot3(d(1),d(2),d(3),'k*')
  plot3(axl(2),d(2),d(3),'k*'), plot3(d(1),axl(4),d(3),'k*'),plot3(d(1),d(2),axl(5),'k*')
  plot3(e(1),e(2),e(3),'go','markersize',radius*1e3), plot3(e(1),e(2),e(3),'g*')
  plot3(axl(2),e(2),e(3),'g*'), plot3(e(1),axl(4),e(3),'g*'),plot3(e(1),e(2),axl(5),'g*')
  plot3(axl(2)*ones(1,tt+1),wmne(i0:i0+tt),wmnt(i0:i0+tt))
  plot3(wmng(i0:i0+tt),axl(4)*ones(1,tt+1),wmnt(i0:i0+tt))
  plot3(wmng(i0:i0+tt),wmne(i0:i0+tt),axl(5)*ones(1,tt+1))
end

%{
figure, plotmatrix([ce(1:end-1);ca(2:end);cae(2:end)]')
figure
plotmatrix([ce(1:end-1);cae(2:end);ca(1:end-1)]')
%}
% plot g, k
figure
fplot(@(x) x+1), hold on, grid on, axis([.46 .47 1.46 1.47])
fplot(@(x) x.^-.5,[0 2])
plot(cg,mnace,'r*')  %plot(.4634,1.469,'r*')

%{
% one-sided A+E pdf
osae = zeros(1,ni);
for i=1:ni
    osae(i) = abs(ce(i))/mnace + ca(randi(ni,1));
end, mean(osae<0)
figure,histogram(osae)

% energy switching sides (v)
mnArn = mean(osae(osae<0))
mnArp = (1-mean(osae(osae<0))*mean(osae<0))/mean(osae>0)
mean(osae(osae<0).^2), mean(osae(osae>0).^2), mean(osae.^2)
mnArn*mean(osae<0)+mnArp*mean(osae>0)
v=-mnArn*mean(osae<0)+mnArp*mean(osae>0)

% approximate kte with gpd
fGpd2 = @(x,th) integral(@(y) fGpd(y,1,2/3,-1/3).*fGpd(x-y,1,2/3,-1/3),-th,th);
FGpd2 = @(x) integral(@(y) fGpd2(y,10),0,x, 'ArrayValued',1)
figure
fplot(@(x) fGpd2(x,10)), grid on, hold on, fplot(@(x) fGpd(x,1,2/3,-1/3)), fplot(@(x) FGpd2(x))

%
figure
winn = 1e4; y = wmne - mean(wmne);
pwelch(y, [], [], [], [], 'twosided')
grid on, set(gca, 'yscale', 'log', 'xscale', 'log')

I=[wmng;wmne;wmnt]';
for i=1:3
    edges = pmi(i):1e-7:pma(i);  %1e-7: 8989,4,3; 1e-8:8989,3,2
    nb=histcounts(I(:,i),edges);
    max(nb)
end

figure
for i=1:ni
    plot3(wmng(1:i),wmne(1:i),wmnt(1:i))
    grid on, axis(axl), hold on
    plot3(wmng(i),wmne(i),wmnt(i), 'r*'), plot3(de,de+1,1,'r*')
    plot3(axl(2)*ones(1,i),wmne(1:i),wmnt(1:i)),plot3(wmng(1:i),axl(4)*ones(1,i),wmnt(1:i)),plot3(wmng(1:i),wmne(1:i),axl(5)*ones(1,i)), 
    hold off
    drawnow %view
    %MV(i) = getframe;    
end% movie(MV)

figure, sif=200;
subplot(131)
pwelch_pl(ce-mean(ce),ni/10,[],[], sif, .95, 'half', 'loglog', 'long-mean')
subplot(132)
pwelch_pl(cae-mean(cae),ni/10,[],[], sif, .95, 'half', 'loglog', 'long-mean')
subplot(133)
pwelch_pl(ca-mean(ca),ni/10,[],[], sif, .95, 'half', 'loglog', 'long-mean')
%set(gca, 'yscale', 'log', 'xscale', 'log')
%}

%% Appendix D: self-consistent transport
% Solving algebraic equations
% symObj = syms; cellfun(@clear,symObj) % delete all syms

% numeric, instead of Cardan's formula, which is cumbersome
roots([1 2 1 -1])
% symbolic
syms x
p0 = x^3 + 2*x^2 + x - 1; 
p1 = v^2*x^3 + 2*v*x^2 + v^2*x - 1;   
p2 = v^2*x^3 + 2*v*x^2 + x - 1; % if v=1.0104, de=.46370315519448729264044650618955
figure,fplot(p0),grid on,hold on,fplot(p1),fplot(p2),axis([.463 .466 -.01 .01])
deroots = solve(p0,x,'MaxDegree', 3); % 'ReturnConditions', true
deroots = simplify(deroots)
de = vpa(deroots(1))
k = 1 + de

% Solving GPD equations
syms sxi ssg sk
S = solve((1-2*sk)*sxi == ssg, ... %sk^2 - 1 == 2*ssg*(1-2*sxi+ssg)/(1-2*sxi)/(1-sxi),...
          (sk^2-1)*(1-2*sxi)*(1-sxi) == 2*ssg*(1-2*sxi+ssg),...
          sxi, ssg, 'ReturnConditions',true);
% solve does not automatically return all solutions of an equation. To return all solutions along with the parameters in the solution and the conditions on the solution, set the ReturnConditions option to true
S2 = solve((1-2*sk)*sxi == ssg, (sk^2-1)*(1-2*sxi)*(1-sxi) == 2*ssg*(1-2*sxi+ssg), ssg, sk, 'ReturnConditions',true);
xi_k = simplify(S.sxi), xi_k = xi_k(4:5)
sg_k = simplify(S.ssg), sg_k = sg_k(4:5)
sg_xi= simplify(S2.ssg)
xi2 = vpa(subs(xi_k,sk,k))
sg2 = vpa(subs(sg_k,sk,k))

figure
subplot(221), ezplot(xi_k(1)), hold on, line([k k],[-2 2],'LineStyle','--')
ezplot(xi_k(2)), grid on, axis([0 2 -2 2]), ylabel('\xi')
subplot(223), ezplot(sg_k(1)), hold on, line([k k],[-1 2],'LineStyle','--')
ezplot(sg_k(2)), grid on, axis([0 2 -1 2]), ylabel('\sigma')
subplot(224), ezplot(sg_xi(1)), hold on
ezplot(sg_xi(2)), ezplot(sg_xi(3)), grid on, axis([-2 2 -1 2])
line([[xi2,xi2]',[-2 -2;2 2]],[-1 -1 sg2'; 2 2 sg2'],'LineStyle','--')

ksi = xi2(2)
sig = sg2(2)
ub = 1 - sig/ksi
ket = 1 + sig/(1-ksi)
kte = k / ket


%% Appendix E: figPfgivenA
setPlotDefaults
SAplusE = @(a,th) (SGpd(th+a,1,2/3,-1/3)+SGpd(th-a,1,2/3,-1/3))/2;

figure
ezplot(@(x) SAplusE(x,1.5))
grid on, axis([-1 1 0 1])
set(gca, 'LineWidth', alw)
xlabel('$a_l$', 'Interpreter', 'latex'), title('')
ylabel('$P(\mathcal{E}_l \neq 0 \; |\; a_l)$', 'Interpreter', 'latex')

%
cf = gcf; set(cf, 'PaperUnits','normalized', 'PaperPositionMode','auto')
print(gcf, '-depsc', '-tiff', 'img/figPfgivenA')

