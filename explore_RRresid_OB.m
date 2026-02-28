function out = explore_RRresid_OB(chanDat, thr)
% explore_RRresid_OB
% Quick-and-dirty exploratory analysis: what OB/gamma features relate to RR_resid?
%
% Usage:
%   S = load(file,'chanDat'); chanDat = S.chanDat;
%   out = explore_RRresid_OB(chanDat, 0.05);
%
% Assumes (per your pipeline):
%   chanDat.use
%   chanDat.behDat.RR_resid
%   chanDat.behDat.length
%   chanDat.gammaBurst (t0_idx_full, t0_sec, prominence, snr, width_sec, phaseID, phaseDist_samp, freqHz)
%   chanDat.gammaBurstSecondary (prominence, snr, width_sec)
%   chanDat.fooof.spectra_flat_log10
%   chanDat.tf.frex
%   chanDat.fooof.gamma_peaks and/or chanDat.fooof.gamma_peak_freq
%   chanDat.targIDX (breaths x 50)
%   chanDat.trial.lowRsp (breaths x time)
%
% No heavy QC/error-checking by design.

if nargin < 2 || isempty(thr), thr = 0.05; end

% ------------------- core vectors -------------------
useVec = chanDat.use(:)==1;

y = double(chanDat.behDat.RR_resid(:));
len = double(chanDat.behDat.length(:));

% Align lengths defensively (minimal)
n = min([numel(useVec), numel(y), numel(len)]);
useVec = useVec(1:n);
y      = y(1:n);
len    = len(1:n);

% ------------------- merge gamma peak matrices (breath x epoch) -------------------
Gpk = merge_gamma_peaks(chanDat);   % breath x 5 (Hz), NaN if none
Gpk = double(Gpk);
Gpk = Gpk(1:n,:);

% ------------------- spectral prominence at the chosen peak (log10 units) -------------------
frex = double(chanDat.tf.frex(:));
flat_log10 = double(chanDat.fooof.spectra_flat_log10);   % breaths x 5 x frex
flat_log10 = flat_log10(1:n,:,:);

gammaBandHz = [25 60];
baselineBandHz = [25 58];
excludeHzAroundPeak = 5;

[promSpec, pkHzSnap, fiPk] = spectral_prominence_from_flat( ...
    flat_log10, frex, Gpk, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak);

% Use epoch 2 as a convenient “canonical” spectral summary (inhale fall)
promSpec_e2 = promSpec(:,2);
pkHz_e2     = pkHzSnap(:,2);

% Also pick each breath’s gammaBurst phaseID epoch if available (1..5)
if isfield(chanDat,'gammaBurst') && isfield(chanDat.gammaBurst,'phaseID')
    phID = double(chanDat.gammaBurst.phaseID(:)); phID = phID(1:n);
else
    phID = nan(n,1);
end

promSpec_byBurstEpoch = nan(n,1);
pkHz_byBurstEpoch     = nan(n,1);
for b = 1:n
    if ~useVec(b), continue; end
    e = phID(b);
    if isfinite(e) && e>=1 && e<=5
        promSpec_byBurstEpoch(b) = promSpec(b,e);
        pkHz_byBurstEpoch(b)     = pkHzSnap(b,e);
    end
end

% ------------------- gammaBurst features -------------------
gb = chanDat.gammaBurst;
t0_sec   = double(gb.t0_sec(:));         t0_sec   = t0_sec(1:n);
t0_full  = double(gb.t0_idx_full(:));    t0_full  = t0_full(1:n);
gb_prom  = double(gb.prominence(:));     gb_prom  = gb_prom(1:n);
gb_snr   = double(gb.snr(:));            gb_snr   = gb_snr(1:n);
gb_wsec  = double(gb.width_sec(:));      gb_wsec  = gb_wsec(1:n);
gb_freq  = double(gb.freqHz(:));         gb_freq  = gb_freq(1:n);

if isfield(gb,'phaseDist_samp')
    gb_phDist = double(gb.phaseDist_samp(:)); gb_phDist = gb_phDist(1:n);
else
    gb_phDist = nan(n,1);
end

% secondary burst “quality” comparison
if isfield(chanDat,'gammaBurstSecondary')
    gs = chanDat.gammaBurstSecondary;
    gs_prom = double(gs.prominence(:));  gs_prom = gs_prom(1:n);
    gs_snr  = double(gs.snr(:));         gs_snr  = gs_snr(1:n);
else
    gs_prom = nan(n,1);
    gs_snr  = nan(n,1);
end

% ------------------- phase at burst: hilbert(lowRsp) AND nearest targIDX (1..50) -------------------
phiHilb = nan(n,1);
phiIdx50 = nan(n,1);
idx50 = nan(n,1);

lowRsp = double(chanDat.trial.lowRsp);          % breaths x time
lowRsp = lowRsp(1:n,:);

targIDX = double(chanDat.targIDX);              % breaths x 50 (full-rate sample indices into trial columns)
targIDX = targIDX(1:n,:);

for b = 1:n
    if ~useVec(b), continue; end
    ii = t0_full(b);
    if ~isfinite(ii), continue; end
    ii = round(ii);
    if ii < 1 || ii > size(lowRsp,2), continue; end

    % hilbert phase (toolbox-free)
    ph = angle(hilbert_fft(lowRsp(b,:).'));
    phiHilb(b) = ph(ii);

    % nearest targIDX (1..50) -> angle
    d = abs(targIDX(b,:) - ii);
    d(~isfinite(d)) = Inf;
    [~,k] = min(d);
    if isfinite(d(k))
        idx50(b) = k;
        phiIdx50(b) = 2*pi*(k-1)/50;  % 0..2pi
    end
end

% ------------------- gamma amplitude at burst (bandpass FFT, hilbert FFT) -------------------
ampGamma = nan(n,1);
for b = 1:n
    if ~useVec(b), continue; end
    ii = t0_full(b);
    if ~isfinite(ii), continue; end
    ii = round(ii);
    if ii < 1 || ii > size(lowRsp,2), continue; end
    f0 = gb_freq(b);
    if ~isfinite(f0), continue; end

    x = double(chanDat.trial.data(b,:)).';
    fs = double(chanDat.fs);
    xb = bandpass_fft(x, fs, max(1,f0-5), min(fs/2-1,f0+5));
    a  = abs(hilbert_fft(xb));
    ampGamma(b) = a(ii);
end

% ------------------- a few “synthetic” predictors -------------------
% a) early vs late burst relative to inhale (t0_sec)
% b) burst “quality” ratio (secondary vs primary)
sec_over_prim = gs_prom ./ gb_prom;

% c) simple “phase consistency” proxy: closeness to phase centroid (lower = tighter)
phCentTight = -gb_phDist;  % higher means closer

% ------------------- predictor table -------------------
T = table();
T.RR_resid = y;
T.breathLen = len;

T.t0_sec = t0_sec;
T.phaseHilb = phiHilb;
T.phaseIdx50 = phiIdx50;
T.idx50 = idx50;

T.gb_prom = gb_prom;
T.gb_snr  = gb_snr;
T.gb_width_sec = gb_wsec;
T.gb_freqHz = gb_freq;

T.sec_prom = gs_prom;
T.sec_snr  = gs_snr;
T.sec_over_prim = sec_over_prim;

T.promSpec_e2 = promSpec_e2;
T.pkHz_e2 = pkHz_e2;
T.promSpec_byBurstEpoch = promSpec_byBurstEpoch;
T.pkHz_byBurstEpoch = pkHz_byBurstEpoch;

T.phID = phID;
T.phCentTight = phCentTight;
T.use = useVec;

% ------------------- quick figures -------------------
idxHi = useVec & isfinite(y) & (y > thr);
idxLo = useVec & isfinite(y) & (y <= thr);

epochNames = ["inhale rise","inhale fall","exhale rise","exhale fall","pause"];

figure('Color','w','Units','normalized','Position',[0.05 0.06 0.92 0.84]);

t = tiledlayout(2,4,'Padding','compact','TileSpacing','compact');
title(t, sprintf('RR_resid exploration | %s | %s | thr=%.3f', string(chanDat.subID), 'chan', thr), ...
    'Interpreter','none','FontWeight','bold');

% (1) your vibe: timing hist high vs low
nexttile(t,1);
histogram(t0_sec(idxHi), -0.5:0.05:2.0, 'Normalization','probability'); hold on;
histogram(t0_sec(idxLo), -0.5:0.05:2.0, 'Normalization','probability');
xlabel('Gamma burst time (s, rel inhale)'); ylabel('Probability');
title('Burst timing: RR\_resid high vs low'); grid on; box off;
legend({sprintf('RR>%.2f (n=%d)',thr,sum(idxHi)), sprintf('RR<=%.2f (n=%d)',thr,sum(idxLo))}, 'Location','best');

% (2) burst phase (hilbert) – quick polar overlay
nexttile(t,2);
polarhistogram(phiHilb(idxHi), 18, 'Normalization','probability'); hold on;
polarhistogram(phiHilb(idxLo), 18, 'Normalization','probability');
title('Hilbert phase at burst'); legend({'high','low'},'Location','southoutside');

% (3) phase via targIDX 1..50 mapping
nexttile(t,3);
polarhistogram(phiIdx50(idxHi), 20, 'Normalization','probability'); hold on;
polarhistogram(phiIdx50(idxLo), 20, 'Normalization','probability');
title('Nearest targIDX phase (1..50 → angle)'); legend({'high','low'},'Location','southoutside');

% (4) RR_resid vs burst timing
nexttile(t,4);
scatter(t0_sec(useVec), y(useVec), 14, 'filled', 'MarkerFaceAlpha',0.45);
xlabel('t0\_sec'); ylabel('RR\_resid');
title('RR\_resid vs burst timing'); grid on; box off;

% (5) RR_resid vs time-domain prominence
nexttile(t,5);
scatter(gb_prom(useVec), y(useVec), 14, 'filled', 'MarkerFaceAlpha',0.45);
xlabel('gammaBurst prominence (time-domain)'); ylabel('RR\_resid');
title('RR\_resid vs burst prominence'); grid on; box off;

% (6) RR_resid vs gamma amplitude at burst (bandpass ±5 Hz)
nexttile(t,6);
scatter(ampGamma(useVec), y(useVec), 14, 'filled', 'MarkerFaceAlpha',0.45);
xlabel('gamma amp at burst (Hilbert)'); ylabel('RR\_resid');
title('RR\_resid vs gamma amplitude'); grid on; box off;

% (7) RR_resid by burst epoch (phaseID)
nexttile(t,7);
if any(isfinite(phID(useVec)))
    boxplot(y(useVec & isfinite(phID)), phID(useVec & isfinite(phID)), 'Labels', cellstr(epochNames));
    ylabel('RR\_resid'); title('RR\_resid by burst epoch'); grid on; box off;
else
    axis off; text(0.1,0.5,'gammaBurst.phaseID missing/NaN', 'Units','normalized');
end

% (8) quick “which features correlate most?” bar (Spearman)
nexttile(t,8);
featNames = ["breathLen","t0_sec","gb_prom","gb_snr","gb_width_sec","gb_freqHz", ...
             "ampGamma","promSpec_byBurstEpoch","sec_over_prim","phCentTight"];
featMat = [T.breathLen, T.t0_sec, T.gb_prom, T.gb_snr, T.gb_width_sec, T.gb_freqHz, ...
           ampGamma, T.promSpec_byBurstEpoch, T.sec_over_prim, T.phCentTight];

rho = nan(1,size(featMat,2));
for k = 1:size(featMat,2)
    rho(k) = spearman_nan(featMat(useVec,k), y(useVec));
end
bar(rho);
set(gca,'XTick',1:numel(featNames),'XTickLabel',featNames,'XTickLabelRotation',30);
ylabel('Spearman \rho'); title('Feature association with RR\_resid');
grid on; box off;

% ------------------- simple multivariate regression (z-scored predictors) -------------------
% Use a modest set to avoid overfitting; you can expand this.
X = [ ...
    zscore_nan(T.breathLen), ...
    zscore_nan(T.t0_sec), ...
    zscore_nan(T.gb_prom), ...
    zscore_nan(T.gb_snr), ...
    zscore_nan(ampGamma), ...
    zscore_nan(T.promSpec_byBurstEpoch), ...
    zscore_nan(T.sec_over_prim) ...
    ];

m = useVec & all(isfinite(X),2) & isfinite(y);

beta = X(m,:) \ y(m);
yhat = X(m,:) * beta;

out = struct();
out.T = T;
out.beta = beta;
out.R2_inSample = 1 - sum((y(m)-yhat).^2) / sum((y(m)-mean(y(m))).^2);

figure('Color','w');
scatter(y(m), yhat, 18, 'filled', 'MarkerFaceAlpha',0.5); grid on; box off;
xlabel('RR\_resid (observed)'); ylabel('RR\_resid (predicted)');
title(sprintf('Linear model (z-scored X) | R^2=%.3f | n=%d', out.R2_inSample, sum(m)));

end

% ======================= helpers =======================

function G = merge_gamma_peaks(chanDat)
% breath x 5, choose gamma_peaks when present; else fall back to gamma_peak_freq per cell
G1 = nan;
G2 = nan;
if isfield(chanDat,'fooof') && isfield(chanDat.fooof,'gamma_peaks') && ~isempty(chanDat.fooof.gamma_peaks)
    G1 = double(chanDat.fooof.gamma_peaks);
end
if isfield(chanDat,'fooof') && isfield(chanDat.fooof,'gamma_peak_freq') && ~isempty(chanDat.fooof.gamma_peak_freq)
    G2 = double(chanDat.fooof.gamma_peak_freq);
end
if ~ismatrix(G1), G1 = nan(size(G2)); end
if ~ismatrix(G2), G2 = nan(size(G1)); end

% reconcile size (minimal)
if any(size(G1) ~= size(G2))
    if numel(G1)==1, G1 = nan(size(G2)); end
    if numel(G2)==1, G2 = nan(size(G1)); end
end

G = G1;
m = ~isfinite(G);
G(m) = G2(m);
end

function [prom, pkHzSnap, fiPk] = spectral_prominence_from_flat(flat_log10, frex, Gpk, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
[nBreaths, nEpochs, nF] = size(flat_log10);
prom = nan(nBreaths, nEpochs);
fiPk = nan(nBreaths, nEpochs);
pkHzSnap = nan(nBreaths, nEpochs);

baseBandMask = frex>=baselineBandHz(1) & frex<=baselineBandHz(2);
gammaMaskAll = frex>=gammaBandHz(1) & frex<=gammaBandHz(2);

for e = 1:nEpochs
    for b = 1:nBreaths
        if ~useVec(b), continue; end
        f0 = Gpk(b,e);
        if ~isfinite(f0), continue; end

        [~,ii] = min(abs(frex - f0));
        fiPk(b,e) = ii;
        pkHzSnap(b,e) = frex(ii);

        y = squeeze(flat_log10(b,e,:));
        if ~any(isfinite(y(gammaMaskAll))), continue; end

        exclMask = frex >= (pkHzSnap(b,e)-excludeHzAroundPeak) & frex <= (pkHzSnap(b,e)+excludeHzAroundPeak);
        bmask = baseBandMask & ~exclMask;

        base = median(y(bmask), 'omitnan');
        if ~isfinite(base), continue; end

        prom(b,e) = y(ii) - base;
    end
end
end

function z = zscore_nan(x)
x = double(x(:));
mu = mean(x(isfinite(x)),'omitnan');
sd = std(x(isfinite(x)),0,'omitnan');
z = (x - mu) ./ sd;
end

function rho = spearman_nan(x,y)
x = double(x(:)); y = double(y(:));
m = isfinite(x) & isfinite(y);
x = x(m); y = y(m);
if numel(x) < 5
    rho = NaN; return
end
rx = tiedrank_simple(x);
ry = tiedrank_simple(y);
rho = corr(rx, ry);
end

function r = tiedrank_simple(x)
% average ranks for ties (toolbox-free)
[xx,ord] = sort(x);
r = zeros(size(x));
n = numel(x);
i = 1;
while i <= n
    j = i;
    while j < n && xx(j+1)==xx(i)
        j = j + 1;
    end
    rk = (i + j) / 2;
    r(ord(i:j)) = rk;
    i = j + 1;
end
end

function xbp = bandpass_fft(x, fs, f1, f2)
% FFT “brick-wall” bandpass, toolbox-free
x = double(x(:));
n = numel(x);
X = fft(x);
f = (0:n-1)' * (fs/n);

mask = (f>=f1 & f<=f2) | (f>=(fs-f2) & f<=(fs-f1));
X(~mask) = 0;
xbp = real(ifft(X));
end

function z = hilbert_fft(x)
% Analytic signal via FFT (toolbox-free hilbert)
x = double(x(:));
n = numel(x);
X = fft(x);
h = zeros(n,1);
if mod(n,2)==0
    h(1) = 1; h(n/2+1) = 1;
    h(2:n/2) = 2;
else
    h(1) = 1;
    h(2:(n+1)/2) = 2;
end
z = ifft(X .* h);
end
