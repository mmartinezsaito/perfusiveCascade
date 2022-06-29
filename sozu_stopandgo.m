% the cascade of Brownian neurons is justified by the 
%  preselection of the generative model
% it's SOC because activities build-up in each layer until firing
% neurons prefer not to fire when there is no signal (FA, hence exn)
% but are forced to not miss to prevent being destroyed (con)
clear, close all
addpath ../Fractals_ScaleInvariance_Patterns/Power-lawDensities_EmpiricalData_09_ClausetShaliziNewman_code/
addpath ../Fractals_ScaleInvariance_Patterns/Ihlen_BeyondScalingLaws/

% define the system
sif = 200; % sampling / input frequency in Hz
prior = 0; 
m = 0; s = 1;
st0 = 3;
ni = 1e6

w = 0.01;
nl = 7

% threshold updating
iseamrand = 0;
eam = 0      % error activity modulation factor:  0 .1 rand   
ueam = eam     % top layer activity suppression factor
lossfunc = 'mmse'  % mmse mmae both

% threshold expansion factor:  0 0.1 
if ~strcmp(lossfunc,'mmse')
    fGpd = @(x, mu, sigma, xi) heaviside(x-mu) .* (1 + xi*(x - mu)/sigma).^(-(xi+1)/xi) ./ sigma .* heaviside(3-x); 
    % average increment for GPD signal
    avgcon = integral(@(x) fGpd(x,1,2/3,-1/3).*(x-1.5), 1.5, 3)/integral(@(x) fGpd(x,1,2/3,-1/3), 1.5, 3);
    tex = avgcon * w % w*3/8: avg quadratic (counterfactual) expansion of theta
    %tex = 0.7 * w      % w*0.7: avg quadratic contraction of theta
    %tex = 0.5 * w      % w*0.5
else
    tex = w         
end
% threshold contraction factor 
tco = w / 1   

% activity suppression rules 
isofw = 0  % is the cascade only forward (or with perfect td suppression)
% bfwud: if after suppression there are still units exceeding threshold, move up, else, down
ismemless = 0   % are activities reset every time step?


% utils
cauchyrnd = @(l, s) l + s * tan(pi * (rand - 1/2));                     % ~1.8
%las1000 = makedist('Stable','alpha',1, 'beta',0, 'gam',1, 'delta',0);  
las0500 = makedist('Stable','alpha',0.5, 'beta',0, 'gam',1, 'delta',0); % ~1.42
% GPD: mu, sigma, ksi
% Pareto: (ksi>0, mu=theta=sigma/ksi) xm=sigma/ksi, alpha=1/ksi 
gp15 = makedist('GeneralizedPareto', 'k',2, 'sigma',2, 'theta',0);    %1/2 * x^(-3/2); ~1.4
gp17 = makedist('GeneralizedPareto', 'k',3/2, 'sigma',3/2, 'theta',0);%2/3 * x^(-5/3); ~1.52
gp20 = makedist('GeneralizedPareto', 'k',1, 'sigma',1, 'theta',0);    %  1 * x^(-2);   ~1.79
gp25 = makedist('GeneralizedPareto', 'k',2/3, 'sigma',2/3, 'theta',0);%3/2 * x^(-5/2); ~2
gp30 = makedist('GeneralizedPareto', 'k',.5, 'sigma',.5, 'theta',0);  %  2 * x^(-3);   ~2.2
gpexp = makedist('GeneralizedPareto', 'k', 0, 'sigma',1, 'theta',0);   %      exp(x)
gpp20 = makedist('GeneralizedPareto', 'k',-1/3, 'sigma',1, 'theta',0); %      x^2
gpp10 = makedist('GeneralizedPareto', 'k',-1/2, 'sigma',1, 'theta',0); %      x
gpp05 = makedist('GeneralizedPareto', 'k',-2/3, 'sigma',1, 'theta',0); %      x^(1/2)
% fractional Brownian motion (Abry & Sellan)
% NS=6 reconstruction steps; sufficiently regular orthogonal wavelet W='db10'
fbm25 = wfbm(0.25, ni);
fbm50 = wfbm(0.5, ni);  % same as Brownian motion
fbm80 = wfbm(0.80, ni);
% Pink noise
pn = pinknoise_hz(1, ni);
% fractional Gaussian noise
fgn05 = diff(wfbm(0.05, ni+1));
fgn25 = diff(wfbm(0.25, ni+1));
fgn50 = diff(wfbm(0.5, ni+1));  % same as randn
fgn80 = diff(wfbm(0.80, ni+1));
fgn95 = diff(wfbm(0.95, ni+1));
% sum of multivariate gaussians
% means and covariance matrices can be simply added
% the expected magnitude of a k-variate gaussian is a k-chi distribution 
kvg_mn = @(k) sqrt(2)*gamma((k+1)/2)/gamma(k/2);
kvg_md = @(k) sqrt(k*(1-2/(9*k))^(1/3));
kvg_ds = @(k) sqrt(k); 
inputpdf = 'fgn50rw'  % 'grw' 'gaussian'      

%%
% initialize variables
A = zeros(ni, nl); a = zeros(1, nl);
T = st0 * ones(ni, nl); st = st0 * ones(1, nl);
E = zeros(ni, nl); 
ACT = zeros(ni, nl); L = [];  
Ebf = []; Abf = [];
i = 1; act = 1; 

while i < ni, i = i + 1; if mod(i, 1e5)==0, i, end
    es = 0; El = zeros(1, nl); 
    
    l = 1; L = [L l];
    
    switch inputpdf
        case 'hgrw',     a(1) = a(1) + abs(normrnd(m, s)); %random('hn', m, s * act);  % Half-normal
        case 'grw',      a(1) = a(1) + normrnd(m, s);  % Gaussian random walk or Wiener process
        case 'cauchyrw', a(1) = a(1) + cauchyrnd(m, s);  % Cauchy process
        case 'gp30rw',   a(1) = a(1) + gp30.random *(2*randi(2)-3);  % GPD (x^-2)
        case 'gpexprw',  a(1) = a(1) + gpexp.random *(2*randi(2)-3);  % GPD (exp(x))
        case 'gpp20rw',  a(1) = a(1) + gpp20.random *(2*randi(2)-3);  % GPD (x^2)
        case 'gpp10rw',  a(1) = a(1) + gpp10.random *(2*randi(2)-3);  % GPD (x)
        case 'unifrw',   a(1) = a(1) + (2*rand-1)*sqrt(3);  
        case 'fgn05rw',  a(1) = a(1) + fgn05(i);
        case 'fgn25rw',  a(1) = a(1) + fgn25(i);
        case 'fgn80rw',  a(1) = a(1) + fgn80(i);
        case 'fgn50rw',  a(1) = a(1) + fgn50(i);
        case 'fgn95rw',  a(1) = a(1) + fgn95(i);
        case 'pnrw',     a(1) = a(1) + pn(i);            
        case 'fbm25rw',  a(1) = a(1) + fbm25(i);
        case 'gaussian', a(1) = normrnd(m, s);   % Gaussian process
        case 'cauchy',   a(1) = cauchyrnd(m, s); % Cauchy process
        case 'las',      a(1) = las0500.random;  % Stable (a=1/2,b=0,l=0,c=1)
        case 'gp',       a(1) = (gp30.random)*(2*randi(2)-3);  % GPD (x^-2)
        otherwise,       break  % don't run the loop if there is no sensory input 
    end
    
    while true
        if iseamrand
            eam = rand; ueam = rand;
        end
        
        % Surprisal random walk / Levy flight
          %a(l) = es;     % random process        
        a(l) = a(l) + es;   % left over A stays 
                        
        if abs(a(l)) <= st(l)   % REST
            if ~isofw && eam~=0, Ebf = [Ebf; zeros(1, nl)]; end
                        
            stl_new = st(l) * (1-tco) + tco * abs(a(l));  % update only after spike volley ends
                                
            % Firing effect (retrospective): error suppression            
            if l > 1
              % update only after spike volley ends
              switch lossfunc
                case 'mmse', stsubl_new = st(1:l-1) * (1-tex) + tex * abs(a(1:l-1)); 
                otherwise,   stsubl_new = min(st(1:l-1) * (1 + tex), abs(a(1:l-1))); 
              end
       
              % top-down effects from updated supraordinate states
              a(1:l-1) = a(1:l-1) * eam;                                      
            end
            if ~isofw && eam~=0, Abf = [Abf; a]; end
                
            % stop spike volley, else continue subvolley from upper or first level 
            if isofw || all(abs(a) <= st) % all(abs(a(1:l)) <= st(1:l)) 
  
                % Resting effect: selective pressure
                % bias increases by energy drift and predictability decreases
                st(l) = stl_new;
    
                if l > 1                     
                    % plasticity/neuromodulation increase predictability, save energy
                    st(1:l-1) = stsubl_new; % update only after spike volley ends
                end
                
                break
            elseif l < nl && any(abs(a(l+1:end)) > st(l+1:end)) 
                % continue from upper level
                l = l + 1;
                es = 0;
            else %if any(abs(a(1:l-1)) > st(1:l-1))  % here, necessarily l>1  
                 % continue from first level (bfw)
                 %l = 1; %l = find(abs(a) > st, 1, 'last');  
                % continue from subordinate level (bfwud)
                l = l - 1;
                es = 0;
            end                   
        else                    % SPIKE
            % error signal
            es = a(l);  
            El(l) = El(l) + es; % abs(es), es         <------SET  
            if ~isofw && eam~=0, Ebf = [Ebf; zeros(1, nl)]; Ebf(end, l) = es; end
            
            if l < nl                
                l = l + 1;     
                if ~isofw && eam~=0, Abf = [Abf; a]; end
            elseif l == nl      
                
                % Firing effects (retrospective): error suppression                 
                switch lossfunc  % update only after spike volley ends 
                  case 'mmse', st1nl_new = st(1:l-1) * (1-tex) + tex * abs(a(1:l-1)); 
                  otherwise,   st1nl_new = min(st(1:l-1) * (1 + tex), abs(a(1:l-1))); 
                end
                % top-down effects from updated supraordinate states                    
                a(1:l-1) = a(1:l-1) * eam;       
   
                % Top layer firing effects (concurrent)
                %  it's unclear what to do in this section
                switch lossfunc  % update only after spike volley ends 
                  case 'mmse', stnl_new = st(l) * (1-tex) + tex * abs(a(l));  
                  otherwise,   stnl_new = min(st(l) * (1 + tex), abs(a(l))); 
                end
                % error suppression 
                a(l) = a(l) * ueam + (1-ueam) * sign(a(l)) * prior; % enable damped itinerancy                                
                % error driving action to sample less randomness 
                %act = st(l) / abs(a(l));   
                if ~isofw && eam~=0, Abf = [Abf; a]; end
                       
                % stop EM loop, else continue subvolley from first level 
                if isofw || all(abs(a) <= st)
                    % plasticity/neuromodulation increase predictability, save energy
                    st(1:nl-1) = st1nl_new;  
  
                    st(nl) = stnl_new;  
                    
                    break
                else %if any(abs(a(1:l-1)) > st(1:l-1)) 
                    % continue from first level (bfw)
                    %l = 1;   %l = find(abs(a) > st, 1, 'last');
                    % continue from subordinate level (bfwud)
                    l = nl - 1;
                    es = 0;
                end                
            end            
        end    
        L = [L l];       
    end    
    
    A(i,:) = a;
    T(i,:) = st;     
    ACT(i) = act;
    E(i,:) = El;
    
   % reset A every time step?
   if ismemless, a = zeros(1, nl); end   
end


%%

figure
[~,edges] = histcounts(log10(T));
[mn, sem, sd] = grpstats(A, [], {'mean','sem', 'std'});
mdT = median(T);
subplot(331), plot(T)
subplot(332), histogram(T, 10.^edges), title(num2str(mdT)), set(gca, 'xscale', 'log')
subplot(333), %plot(ACT)
subplot(334), plot(A), hold on, plot(sum(A,2))
subplot(335), errorbar(mn, sem), hold on, plot(sd)
subplot(336), histogram(sum(abs(A),2)), hold on, histogram(sum(A,2))
subplot(3,3,7:8), plot(Abf), hold on, plot(sum(Abf,2)), plot(L'*10, '.:') 
subplot(339), histogram(sum(abs(Abf),2)), hold on, histogram(sum(Abf,2))
ss = sprintf('img/bfwud_input%s_lf%s_nl%02u_ni%u_w%.2f_isrand%u_eam%.1f_ueam%.1f', inputpdf, lossfunc, nl, ni, w, iseamrand, eam, ueam);
preprint


%%

% Spectral analysis
% compare raw and abs spectra
%{
figure, S = L; Sbf = Lbf;
subplot(221), for i=1:nl, pwelch_pl(S(:,i), [], [], [], sif), hold on, end, grid on
for i=1:nl, pwelch_pl(abs(S(:,i)), [], [], [], sif), end,  set(gca, 'yscale', 'log', 'xscale', 'log')
subplot(222), pwelch_pl(sum(S,2), [], [], [], sif), grid on, hold on
pwelch_pl(sum(abs(S),2), [], [], [], sif), set(gca, 'yscale', 'log', 'xscale', 'log')
subplot(223), for i=1:nl, pwelch_pl(Sbf(:,i), [], [], [], sif), hold on, end, grid on
for i=1:nl, pwelch_pl(abs(Sbf(:,i)), [], [], [], sif), end,  set(gca, 'yscale', 'log', 'xscale', 'log')
subplot(224), pwelch_pl(sum(Sbf,2), [], [], [], sif), grid on, hold on
pwelch_pl(sum(abs(Sbf),2), [], [], [], sif), set(gca, 'yscale', 'log', 'xscale', 'log')

pwelch([Ebf sum(Ebf,2)], [], [], [], sif, 'ConfidenceLevel', .95, 'psd', 'onesided') ,%set(gca, 'xscale', 'log')
% pcov pmcov pyulear pburg peig pmusic pmtm
%}


% scaling exponent DFA alpha for topmost layer, g=2-2a, b=2a-1, g=1-b
for i=1:nl
    %thresholds
    a = DFA_fun(T(:,i), [10 ni/10]);  % a=1.5 (b=2)
    H(i,1) = a(1); 
    % top-down activity
    a = DFA_fun(A(:,i), [10 ni/10]);  % a=1.34 (b=1.7)
    H(i,2) = a(1); 
    % bottom-up activity
    a = DFA_fun(E(:,i), [10 ni/10]);  % a= (b=)
    H(i,3) = a(1); 
    a = DFA_fun(sum(E(:,1:i),2), [10 ni/10]);  % a= (b=)
    H(i,4) = a(1); 
    %td+bu
    a = DFA_fun(sum([A(:,1:i) E(:,1:i)],2), [10 ni/10]); H(i,5)=a(1); 
end
a = DFA_fun(sum(T,2), [10 ni/10]); H(nl+1,1) = a(1); 
a = DFA_fun(sum(A,2), [10 ni/10]); H(nl+1,2) = a(1); 
a = DFA_fun(sum(E,2), [10 ni/10]); H(nl+1,3) = a(1); 
a = DFA_fun(sum([A E],2), [10 ni/10]); H(nl+1,5) = a(1); 
H    


%% Multifractal analysis
figure
nr = ceil(sqrt(nl+9));
lc = {'L1','L2','L3','L4','L5','L6','L7','sumL','wgn', 'pn', 'fbmH0', 'bm', 'fbmH1', 'E', 'sumLbf', 'sumEbf'};
for i=1:nl+9
    switch i
        case nl+1, ts = sum(A,2);
        case nl+2, ts = randn(1,ni);
        case nl+3, ts = pinknoise_hz(1,ni);
        case nl+4, ts = wfbm(0,ni);      % Hoder exponent H= 0
        case nl+5, ts = wfbm(.5,ni);     % Hoder exponent H= .5
        case nl+6, ts = wfbm(1,ni);      % Hoder exponent H= 1
        case nl+7, ts = sum(E,2);
        case nl+8, ts = sum(Abf,2);
        case nl+9, ts = sum(Ebf,2);
        otherwise, ts = A(:,i); 
    end
    
    % Power spectral density estimation with Welch's method    
    subplot(nr,nr,i)
    [Pxx,F] = pwelch(ts);
    plot(log10(F(2:end)),log10(Pxx(2:end))), grid on, hold on      
    Xpred = [ones(length(F(2:end)),1) log10(F(2:end))];
    b = lscov(Xpred,log10(Pxx(2:end)));
    y = b(1)+b(2)*log10(F(2:end));
    plot(log10(F(2:end)),y,'r--')
    xlabel('log10(F)'), ylabel('log10(Pxx)')
    title([lc{i} ' PSD slope estimate: ' num2str(b(2))])
    
    % WTMM
    [hexp(i) tauq_wm(i,:) sf_wm(i)] = wtmm(ts, 'ScalingExponent', 'global');
    %[localhexp, wt, wavscales] = wtmm(ts, 'ScalingExponent', 'local'); % for finding local discontinuities
    fprintf('WTMM estimate is -2*(%1.2f)-1 = %1.2f\n',hexp(i),-2*hexp(i)-1);   
    
    % Wavelet leaders are derived from the critically sampled DWT coefficients. 
    %WLs offer significant theoretical advantages over wavelet coefficients in 
    %the multifractal formalism. WL are time- or space-localized suprema of the 
    %absolute value of the DW coefficients. The time localization of the suprema 
    %requires that DW coefficients are obtained using a compactly supported 
    %wavelet. The Holder exponents, which quantify the local regularity, are 
    %determined from these suprema. The singularity spectrum indicates the size 
    %of the set of Holder exponents in the data.
    [dh(i,:) he(i,:) cp(i,:) tauq_wl(i,:) wl sf_wl(i)] = dwtleader(ts);
    % singularity spectrum, Holder exp, log-cumulants, scaling exp, WL, multiresolution sf
    fprintf('Wavelet leader estimate is -2*(%1.2f)-1 = %1.2f\n',cp(i,1),-2*cp(i,1)-1);    
end

% plot multifractal spectrum D(h) vs h
%The singularity spectrum Dh is estimated using structure functions determined 
%for the linearly-spaced moments in [–5:5]. The structure functions are 
%computed based on the wavelet leaders obtained using the biorthogonal spline 
%wavelet filter. The closer a Holder exponent is to 1 (0), the closer the 
%function is to differentiable (discontinuous). 
figure
for i=1:nl+9
    ph = plot(he(i,:), dh(i,:), '-o'); ph.LineWidth = 1;
    switch i
        case 1,                ph.LineStyle = '-.'; hold on
        case {2,3,4,5,6,7},    ph.LineStyle = '--'; ph.Marker = 'none';
        case nl+1,             set(ph, 'LineStyle','-', 'Marker', '+')
        case nl+2,             ph.LineStyle = ':'; ph.Marker = '^';
        case nl+3,             ph.LineStyle = ':'; ph.Marker = '^';
        case {nl+4,nl+5,nl+6}, ph.LineStyle = ':'; ph.Marker = 'p';
        case nl+7,             ph.Marker = '*';
        case nl+8,             ph.LineStyle = '-.'; ph.Marker = '+';
        case nl+9,             ph.LineStyle = '-.'; ph.Marker = '*';
        %otherwise, set(ph, 'LineStyle','-.', 'Marker', '+')        
    end
end
title('Multifractal Spectrum'), xlabel('h'), ylabel('D(h)')
legend(lc), grid on
preprint


% A monofractal process has a linear scaling law as a function of the 
%statistical moments, while a multifractal process has a nonlinear scaling 
%law. dwtleader uses the range of moments from -5 to 5 in estimating these 
%scaling laws. Linearity is indicated by zeroish 2nd and 3rd cumulants. The 
%1st cumulant is the slope of scaling exponents vs moments
figure
plot(-5:5, tauq_wl(nl,:), 'bo--'), title('Estimated scaling exponents');
grid on, xlabel('qth Moments'), ylabel('\tau(q)')
fprintf('Cumulants: %1.2f %1.2f %1.2f\n', cp(nl,:))

%
% up to layer i
for i=1:nl+2
    switch i
        case nl+1, tsa = sum(abs(Ebf), 2); 
        case nl+2, tsa = avlen; 
        otherwise, tsa = sum(abs(E(:,1:i)),2); 
    end
    [esia(i), lb(i), ll] = plfit(tsa(tsa~=0), 'limit', 10, 'nosmall');     
    if i >= nl, plplot(tsa(tsa~=0), lb(i), esia(i)); grid on, end
    %preprint
end, esia
%plvar(tsa)
%plpva(tsa, lb)
%[besia, blb, bll] = plfit(B); plplot(B, blb, besia)


%%
% plot (partition function) scaling exponents as a function of moments
ph = plot(-2:0.1:2, tauq_wm(i,:), 'b-o'); grid on
title('Monofractal (or multifractal) scaling exponents')
xlabel('Qth Moment')
ylabel('Scaling Exponents')
% plot wavelet maxima lines and estimates of local Holder exponents
wtmm(ts, 'ScalingExponent', 'local')
% plot tauq vs. structfunc
figure
x = ones(length(sf(i).logscales),2); x(:,2) = sf(i).logscales;
betahat = lscov(x, sf(i).Tq, sf(i).weights); betahat = betahat(2,:);
subplot(121)
plot(-2:.1:2,tauq(i,:)), grid on
title('From tauq output'), xlabel('Qth Moment'), ylabel('Scaling Exponents')
subplot(122)
plot(-2:.1:2, betahat(1:41)), grid on
title('From structfunc output'), xlabel('Qth Moment')

%%
figure, plotmatrix(A)
figure, plotmatrix(T)
figure, plotmatrix(E)
figure, plotmatrix(A,T)
figure, plotmatrix(A,E)
figure, plotmatrix(T,E)

%
figure
subplot(221), autocorr(sum(A,2), ni/10)
subplot(222), autocorr(ts, ni/10)
subplot(223), crosscorr(A(:,1), A(:,nl), ni/10) % xcorr
subplot(224), autocorr(ACT, ni/10)

figure
subplot(221), autocorr(A(:, 3), 30)
subplot(222), crosscorr(A(:,1), A(:,nl), 30)
subplot(223), parcorr(A(:, 3))
subplot(224), crosscorr(A(:,4), A(:,nl), 30)

subplot(121), periodogram(A(:, [1 nl]))
subplot(122), for i=1:nl, histogram(A(:,i)), hold on, end
