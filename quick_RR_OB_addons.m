% ===================== Add these two plots to your workflow =====================
% 1) 2D binned heatmap: burst phase bin × burst prominence (color = mean RR_resid)
% 2) Median-split RR_resid: ITPC from gammaLockTF.phase_primary (breaths × time × freq)

function [figHeat, figITPC, out] = quick_RR_OB_addons(chanDat, nPhaseBins, nPromBins)
if nargin < 2 || isempty(nPhaseBins), nPhaseBins = 50; end     % matches targIDX resolution
if nargin < 3 || isempty(nPromBins),  nPromBins  = 8;  end     % coarse bins for readability

useVec = chanDat.use(:)==1;
RR = double(chanDat.behDat.RR_resid(:));

n = min([numel(useVec), numel(RR)]);
useVec = useVec(1:n);
RR     = RR(1:n);

% ---------- Burst phase bins (two definitions) ----------
% A) hilbert phase from lowRsp at t0_idx_full
% B) targIDX-nearest bin (1..50) -> angle (0..2pi)
t0_full = double(chanDat.gammaBurst.t0_idx_full(:)); t0_full = t0_full(1:n);

phiHilb = nan(n,1);
phiIdx50 = nan(n,1);
idx50 = nan(n,1);

lowRsp = double(chanDat.trial.lowRsp); lowRsp = lowRsp(1:n,:);
targIDX = double(chanDat.targIDX);     targIDX = targIDX(1:n,:);

for b = 1:n
    if ~useVec(b), continue; end
    ii = round(t0_full(b));
    if ~isfinite(ii) || ii<1 || ii>size(lowRsp,2), continue; end

    ph = angle(hilbert_fft(lowRsp(b,:).'));
    phiHilb(b) = ph(ii);

    d = abs(targIDX(b,:) - ii);
    d(~isfinite(d)) = Inf;
    [~,k] = min(d);
    if isfinite(d(k))
        idx50(b) = k;
        phiIdx50(b) = 2*pi*(k-1)/50;  % 0..2pi
    end
end

% ---------- Burst prominence (time-domain) ----------
prom = double(chanDat.gammaBurst.prominence(:)); prom = prom(1:n);

% ===================== (1) 2D heatmap: phase × prominence =====================
% Use targIDX phase by default (more “respiratory-phase aligned” in your pipeline)
phaseAng = phiIdx50;
figHeat = plot_phase_prom_heatmap(RR, prom, phaseAng, useVec, nPhaseBins, nPromBins);

% ===================== (2) Median split RR_resid -> ITPC from phase_primary =====================
figITPC = plot_ITPC_medianSplit_RR(chanDat, RR, useVec);

% Return a few things for downstream use
out = struct();
out.useVec = useVec;
out.RR = RR;
out.prom = prom;
out.phiHilb = phiHilb;
out.phiIdx50 = phiIdx50;
out.idx50 = idx50;
end

% -------------------------------------------------------------------------
function fig = plot_phase_prom_heatmap(RR, prom, phaseAng, useVec, nPhaseBins, nPromBins)

m = useVec & isfinite(RR) & isfinite(prom) & isfinite(phaseAng);

% Phase edges: cover 0..2pi (since phiIdx50 is 0..2pi). If you pass hilbert phase, it’s -pi..pi;
% still fine—just swap edges below if you prefer.
phMin = 0; phMax = 2*pi;
if any(phaseAng(m) < 0) % likely hilbert phase
    phMin = -pi; phMax = pi;
end
phaseEdges = linspace(phMin, phMax, nPhaseBins+1);

% Prominence edges (quantile bins for stability)
promEdges = quantile_edges_simple(prom(m), nPromBins);

ip = discretize(phaseAng(m), phaseEdges);
ib = discretize(prom(m), promEdges);

ok = isfinite(ip) & isfinite(ib);
ip = ip(ok); ib = ib(ok);
rr = RR(m); rr = rr(ok);

% Accumulate mean RR per (phaseBin, promBin)
lin = sub2ind([nPhaseBins, nPromBins], ip, ib);
sumRR = accumarray(lin, rr, [nPhaseBins*nPromBins 1], @sum, 0);
cnt   = accumarray(lin, 1,  [nPhaseBins*nPromBins 1], @sum, 0);
meanRR = sumRR ./ max(cnt,1);

M = reshape(meanRR, [nPhaseBins, nPromBins]);
C = reshape(cnt,    [nPhaseBins, nPromBins]);

% Centers for axis labels
phaseCtr = 0.5*(phaseEdges(1:end-1) + phaseEdges(2:end));
promCtr  = 0.5*(promEdges(1:end-1)  + promEdges(2:end));

fig = figure('Color','w','Units','normalized','Position',[0.08 0.10 0.86 0.72]);
t = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
title(t,'RR\_resid as function of burst phase × burst prominence','FontWeight','bold');

% Mean RR heatmap
nexttile(t,1);
imagesc(promCtr, phaseCtr, M);
set(gca,'YDir','normal');
xlabel('Burst prominence bin (time-domain)');
ylabel('Burst phase (rad)');
title('Mean RR\_resid in each bin');
colorbar; grid on; box off;

% Counts heatmap
nexttile(t,2);
imagesc(promCtr, phaseCtr, C);
set(gca,'YDir','normal');
xlabel('Burst prominence bin (time-domain)');
ylabel('Burst phase (rad)');
title('Counts per bin');
colorbar; grid on; box off;

end

% -------------------------------------------------------------------------
function fig = plot_ITPC_medianSplit_RR(chanDat, RR, useVec)

m0 = useVec & isfinite(RR);
medRR = median(RR(m0), 'omitnan');

idxHi = m0 & RR >= medRR;
idxLo = m0 & RR <  medRR;

% breaths × time × freq (your stated convention)
ph = double(chanDat.gammaLockTF.phase_primary);

% Try to pull axes; if missing, just use indices
if isfield(chanDat.gammaLockTF,'tVec')
    tVec = double(chanDat.gammaLockTF.tVec(:));
else
    tVec = (1:size(ph,2)).';
end
if isfield(chanDat.gammaLockTF,'frexSel')
    fVec = double(chanDat.gammaLockTF.frexSel(:));
elseif isfield(chanDat.gammaLockTF,'frex')
    fVec = double(chanDat.gammaLockTF.frex(:));
else
    fVec = (1:size(ph,3)).';
end

% If dimensions don’t match (common gotcha), do a tiny “best guess” permute:
% We want: [breaths, time, freq]
if size(ph,2) ~= numel(tVec) && size(ph,3) == numel(tVec)
    ph = permute(ph, [1 3 2]); % breaths, time, freq
end
if size(ph,3) ~= numel(fVec) && size(ph,2) == numel(fVec)
    ph = permute(ph, [1 3 2]); % breaths, time, freq
end

itpc_hi = squeeze(abs(mean(exp(1i*ph(idxHi,:,:)), 1, 'omitnan'))); % time × freq
itpc_lo = squeeze(abs(mean(exp(1i*ph(idxLo,:,:)), 1, 'omitnan'))); % time × freq
ditpc   = itpc_hi - itpc_lo;

% Plot
fig = figure('Color','w','Units','normalized','Position',[0.06 0.08 0.90 0.78]);
t = tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');
title(t, sprintf('ITPC from phase\\_primary | RR median split (med=%.4f) | nHi=%d nLo=%d', ...
    medRR, sum(idxHi), sum(idxLo)), 'FontWeight','bold');

% Hi
nexttile(t,1);
imagesc([], [], itpc_hi'); set(gca,'YDir','normal');
yticks(10:10:100)
yticklabels(round(fVec(10:10:100)))
xlabel('Time'); ylabel('Frequency (Hz)');
title('ITPC | RR\_high'); colorbar; grid on; box off;
clim([0 .2]);

% Lo
nexttile(t,2);
imagesc(tVec, [], itpc_lo'); set(gca,'YDir','normal');
yticks(10:10:100)
yticklabels(round(fVec(10:10:100)))
xlabel('Time'); ylabel('Frequency (Hz)');
title('ITPC | RR\_low'); colorbar; grid on; box off;
clim([0 .2]);

% Diff
nexttile(t,3);
imagesc(tVec, [], ditpc'); set(gca,'YDir','normal');
yticks(10:10:100)
yticklabels(round(fVec(10:10:100)))
xlabel('Time'); ylabel('Frequency (Hz)');
title('ITPC diff (high - low)'); colorbar; grid on; box off;
clim([-0.2 0.2]);

end

% -------------------------------------------------------------------------
function edges = quantile_edges_simple(x, nBins)
% Toolbox-free quantile bin edges (monotone, covers full range)
x = sort(x(:));
x = x(isfinite(x));
if isempty(x)
    edges = linspace(0,1,nBins+1);
    return
end
p = linspace(0,1,nBins+1);
idx = 1 + (numel(x)-1)*p;
lo = floor(idx); hi = ceil(idx);
lo(lo<1)=1; hi(hi<1)=1; lo(lo>numel(x))=numel(x); hi(hi>numel(x))=numel(x);
frac = idx - lo;
edges = x(lo) + frac(:).*(x(hi)-x(lo));

edges(1) = min(x);
edges(end) = max(x);

% enforce strictly increasing edges (discretize hates duplicates)
for i = 2:numel(edges)
    if edges(i) <= edges(i-1)
        edges(i) = edges(i-1) + eps(edges(i-1) + 1);
    end
end
end

% -------------------------------------------------------------------------
function z = hilbert_fft(x)
% Analytic signal via FFT (toolbox-free Hilbert)
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
