function run_gamma_figures_pipeline()
% RUN_GAMMA_FIGURES_PIPELINE
% Participant/session pipeline:
%   1) group chanDat files by stem (filename without final _NN)
%   2) gammaQC across channels, choose winner
%   3) generate 3 figure sets (split into multiple multipanel figures)
%
% Assumptions (per your clarifications):
% - Each processed file is a chanDat struct saved as variable 'chanDat'
% - Filename ends with _NN.mat where NN is 2-digit channel number
% - Processed files live ONLY in the top level of CHANDAT_processed (ignore subdirs)
% - Fields exist: chanDat.use, behDat.length, fooof.gamma_peaks, fooof.gamma_peak_freq,
%   fooof.spectra_flat_log10, fooof.aperiodic_params, tf.frex, gammaBurst, gammaBurstSecondary,
%   gammaLockTF.*
%
% No toolboxes required for KDE; filtering/analytic signal are FFT-based here.

% ------------------- USER PATHS -------------------
procDir = "R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\CHANDAT_processed\";
outRoot = "R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\SignalProcessing\tf_plots\";

% Subdirectories
dirGammaQC        = fullfile(outRoot, "gammaQC");
dirGammaOsEvidence= fullfile(outRoot, "GammaOsEvidence");
dirBurstTiming    = fullfile(outRoot, "BurstTiming");
dirGammaEventStruct=fullfile(outRoot,"GammaEventStruct");

mkDirIfNeeded(dirGammaQC);
mkDirIfNeeded(dirGammaOsEvidence);
mkDirIfNeeded(dirBurstTiming);
mkDirIfNeeded(dirGammaEventStruct);

% ------------------- GLOBAL SETTINGS -------------------
epochNames = ["inhale rise","inhale fall","exhale rise","exhale fall","pause"];
nEpochs = 5;

% Prominence baseline settings (FOOOF flattened log10)
gammaBandHz     = [25 60];
baselineBandHz  = [25 58];
excludeHzAroundPeak = 5;

% Event-locked ERP settings
erpWinSec = 1.5;
bandHalfWidthHz = 5;  % bandpass = f0 ± 5 Hz

% GammaLockTF plotting limits
zPowCLim  = [-5 10];
itpcCLim  = [-0.2 0.2];

% Figure export settings
dpi = 200;

% ------------------- FIND & GROUP FILES -------------------
files = dir(fullfile(procDir, "*.mat"));
files = files(~[files.isdir]);
if isempty(files)
    error("No .mat files found in %s", procDir);
end

% Group key = filename minus final _NN
keys = strings(numel(files),1);
for i = 1:numel(files)
    keys(i) = stripChanSuffix(files(i).name);
end

[uKeys, ~, gidx] = unique(keys, 'stable');

fprintf("Found %d session-groups.\n", numel(uKeys));

% ------------------- MAIN LOOP OVER SESSION GROUPS -------------------
for gi = 1:numel(uKeys)
    groupKey = uKeys(gi);
    idx = find(gidx == gi);
    groupFiles = files(idx);

    fprintf("\n=== Group %d/%d: %s  (%d channels) ===\n", gi, numel(uKeys), groupKey, numel(groupFiles));

    % Load all channels in group and compute gammaQC metrics
   chanSumm = struct([]);  % <-- NOT repmat(struct(),...)

    for ci = 1:numel(groupFiles)
        fpath = fullfile(groupFiles(ci).folder, groupFiles(ci).name);
        S = load(fpath, "chanDat");
        chanDat = S.chanDat;
    
        s = summarize_channel_gamma(chanDat, epochNames, ...
            gammaBandHz, baselineBandHz, excludeHzAroundPeak);
    
        s.filePath = string(fpath);
        s.fileName = string(groupFiles(ci).name);
    
        if ci == 1
            chanSumm = repmat(s, numel(groupFiles), 1);  % template has all fields
        end
    
        chanSumm(ci) = s;
    end

    % Pick winner (robust z-scored composite within this group)
    [winnerIdx, scoreTable, chanSumm] = pick_winner_channel(chanSumm);

    winner = chanSumm(winnerIdx);
    fprintf("Winner: %s  (score=%.3f)\n", winner.chanLabelSafe, winner.score);

    % Build a stable prefix for filenames
    sessNum = winner.sessNum;
    taskName = chanDat.task;
    prefix = sprintf("%s_sess%d_%s_%s", winner.subIDSafe, sessNum, winner.chanLabelSafe, taskName);

    % ------------------- STAGE 1: gammaQC figure -------------------
    fig = plot_gammaQC_group(scoreTable, epochNames, groupKey);
    outFile = fullfile(dirGammaQC, prefix + "_gammaQC.jpg");
    exportgraphics(fig, outFile, "Resolution", dpi);
    close(fig);

    % ------------------- LOAD WINNER chanDat once -------------------
    Sw = load(winner.filePath, "chanDat");
    chanDatW = Sw.chanDat;

    % ------------------- STAGE 2: GammaOsEvidence (split into 2 figs) -------------------
    % gamma oscilation peak frequency and prominence
    
    fig = plot_GammaOsEvidence_A1A2(chanDatW, epochNames, ...
        gammaBandHz, baselineBandHz, excludeHzAroundPeak);
    outFile = fullfile(dirGammaOsEvidence, prefix + "_GammaOsEvidence_A1A2.jpg");
    exportgraphics(fig, outFile, "Resolution", dpi);
    close(fig);
   
    % gamma oscilation peak + prominence across length and condition
    % try
    % fig = plot_GammaOsEvidence_A3A4(chanDatW, epochNames, ...
    %     gammaBandHz, baselineBandHz, excludeHzAroundPeak);
    % outFile = fullfile(dirGammaOsEvidence, prefix + "_GammaOsEvidence_A3A4.jpg");
    % exportgraphics(fig, outFile, "Resolution", dpi);
    % close(fig);
    % catch
    % end
    % ------------------- STAGE 3: BurstTiming (split into 2 figs) -------------------
    % B1 + B2
    
    fig = plot_BurstTiming_B1B2(chanDatW, epochNames);
    outFile = fullfile(dirBurstTiming, prefix + "_BurstTiming_B1B2.jpg");
    exportgraphics(fig, outFile, "Resolution", dpi);
    close(fig);
  
    % 
    % % B3 + B4
    % try
    % fig = plot_BurstTiming_B3B4(chanDatW, epochNames, erpWinSec, bandHalfWidthHz);
    % outFile = fullfile(dirBurstTiming, prefix + "_BurstTiming_B3B4.jpg");
    % exportgraphics(fig, outFile, "Resolution", dpi);
    % close(fig);
    % catch
    % end
    % 
    % % ------------------- STAGE 4: GammaEventStruct (2 figs) -------------------
    % % Overall + length split
    % try
    % fig = plot_GammaEventStruct_Overall_LenSplit(chanDatW, zPowCLim, itpcCLim);
    % outFile = fullfile(dirGammaEventStruct, prefix + "_GammaEventStruct_Overall_LenSplit.jpg");
    % exportgraphics(fig, outFile, "Resolution", dpi);
    % close(fig);
    % catch
    % end
    % 
    % % Quality split
    % fig = plot_GammaEventStruct_QualitySplit(chanDatW, zPowCLim, itpcCLim);
    % outFile = fullfile(dirGammaEventStruct, prefix + "_GammaEventStruct_QualitySplit.jpg");
    % exportgraphics(fig, outFile, "Resolution", dpi);
    % close(fig);

end

fprintf("\nDone.\n");

% ========================= HELPERS =========================

    function mkDirIfNeeded(d)
        if ~exist(d, "dir"), mkdir(d); end
    end

end % main function



% ============================================================
% ===================== CORE SUMMARY =========================
% ============================================================

function summ = summarize_channel_gamma(chanDat, epochNames, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
% Summarize gamma characteristics for channel selection

useVec = chanDat.use(:)==1;
nBreaths = numel(useVec);
nEpochs = numel(epochNames);

summ.subID   = string(chanDat.subID);
summ.subIDSafe = sanitize_for_filename(summ.subID);

if isfield(chanDat,'sessNum') && ~isempty(chanDat.sessNum) && isfinite(chanDat.sessNum)
    summ.sessNum = double(chanDat.sessNum);
else
    summ.sessNum = 1;
end

% Channel label for filenames
try
    cl = string(chanDat.labels{chanDat.chi});
catch
    cl = "chan" + string(chanDat.chi);
end
summ.chanLabel = cl;
summ.chanLabelSafe = sanitize_for_filename(cl);

% ----- gamma peak detection matrix -----
G = chanDat.fooof.gamma_peaks; % breaths x 5 (NaN = no peak)
G = double(G);

% Prevalence (QC breaths only)
prev = nan(1,nEpochs);
cnt  = zeros(1,nEpochs);
den  = sum(useVec);
for e = 1:nEpochs
    good = useVec & isfinite(G(:,e));
    cnt(e) = sum(good);
    prev(e) = cnt(e) / max(den,1);
end
summ.prev_by_epoch = prev;
summ.count_by_epoch = cnt;

% Frequency MAD (QC breaths only, pooled and per epoch)
madEpoch = nan(1,nEpochs);
for e=1:nEpochs
    x = G(useVec & isfinite(G(:,e)), e);
    madEpoch(e) = robust_mad(x);
end
summ.madFreq_by_epoch = madEpoch;

xAll = G(useVec & isfinite(G));
summ.madFreq_pooled = robust_mad(xAll);

% ----- prominence (QC breaths only; pooled + per epoch) -----
[promEpoch, promAll] = compute_gamma_prominence(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak);
summ.prom_by_epoch = promEpoch;
summ.prom_pooled   = promAll;

% ----- burst quality (QC breaths only) -----
if isfield(chanDat,'gammaBurst') && isfield(chanDat.gammaBurst,'prominence')
    gbProm = double(chanDat.gammaBurst.prominence(:));
    gbSNR  = [];
    if isfield(chanDat.gammaBurst,'snr')
        gbSNR = double(chanDat.gammaBurst.snr(:));
    end
    have = useVec & isfinite(gbProm);
    summ.burstProm_median = median(gbProm(have), 'omitnan');
    summ.burstCoverage    = sum(have)/max(sum(useVec),1);
    if ~isempty(gbSNR)
        summ.burstSNR_median = median(gbSNR(useVec & isfinite(gbSNR)), 'omitnan');
    else
        summ.burstSNR_median = NaN;
    end
else
    summ.burstProm_median = NaN;
    summ.burstCoverage    = 0;
    summ.burstSNR_median  = NaN;
end

summ.score = NaN; % filled later
summ.filePath = "";
summ.fileName = "";

end


function [winnerIdx, T, chanSumm] = pick_winner_channel(chanSumm)
% Robust within-group scoring to select winner

n = numel(chanSumm);
chanLabel = strings(n,1);
subID = strings(n,1);
sessNum = nan(n,1);

prevPooled = nan(n,1);
promPooled = nan(n,1);
madFreq    = nan(n,1);
burstProm  = nan(n,1);
burstCov   = nan(n,1);

for i=1:n
    chanLabel(i) = chanSumm(i).chanLabelSafe;
    subID(i)     = chanSumm(i).subIDSafe;
    sessNum(i)   = chanSumm(i).sessNum;

    prevPooled(i) = mean(chanSumm(i).prev_by_epoch, 'omitnan');
    promPooled(i) = chanSumm(i).prom_pooled;
    madFreq(i)    = chanSumm(i).madFreq_pooled;
    burstProm(i)  = chanSumm(i).burstProm_median;
    burstCov(i)   = chanSumm(i).burstCoverage;
end

% Robust z (median/MAD)
zPrev  = robust_z(prevPooled);
zProm  = robust_z(promPooled);
zMadF  = robust_z(madFreq);
zBProm = robust_z(burstProm);
zBCov  = robust_z(burstCov);

score = zPrev + zProm;

% Package table sorted
T = table(chanLabel, subID, sessNum, prevPooled, promPooled, madFreq, burstProm, burstCov, ...
    zPrev, zProm, zMadF, zBProm, zBCov, score);

[~, ord] = sort(score, 'descend', 'MissingPlacement','last');
T = T(ord,:);

% Winner = first row
winnerChanLabel = T.chanLabel(1);

% map back to index in chanSumm
winnerIdx = find(strcmp(chanLabel, winnerChanLabel), 1, 'first');

% Save score into chanSumm (optional)
for i=1:n
    chanSumm(i).score = score(i);
end

end



% ============================================================
% ===================== FIGURES: gammaQC =====================
% ============================================================

function fig = plot_gammaQC_group(T, epochNames, groupKey)
fig = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.85]);
t = tiledlayout(fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');
title(t, "gammaQC: " + string(groupKey), 'Interpreter','none','FontWeight','bold');

% Panel 1: composite score
nexttile(t,1);
bar(T.score);
set(gca,'XTick',1:height(T),'XTickLabel',T.chanLabel,'XTickLabelRotation',45);
ylabel('Composite score (robust z-sum)');
grid on; box off;
title('Channel ranking (higher = better)');
hold on;
plot(1, T.score(1), 'kp', 'MarkerFaceColor','k', 'MarkerSize',10);
hold off;

% Panel 2: prevalence by epoch heatmap
nexttile(t,2);
% Build matrix channels x epochs
nC = height(T);
nE = numel(epochNames);
Mprev = nan(nC,nE);
for i=1:nC
    % Not stored per epoch in T; approximate using pooled only for QC figure.
    % We'll show pooled prevalence bar instead of epoch-heatmap here for simplicity.
end
% Use pooled prevalence bar
bar(T.prevPooled);
set(gca,'XTick',1:nC,'XTickLabel',T.chanLabel,'XTickLabelRotation',45);
ylabel('Proportion (QC breaths)');
grid on; box off;
title('Peak prevalence (pooled across epochs)');

% Panel 3: prominence pooled
nexttile(t,3);
bar(T.promPooled);
set(gca,'XTick',1:nC,'XTickLabel',T.chanLabel,'XTickLabelRotation',45);
ylabel('Median prominence (log10 flat residual)');
grid on; box off;
title('Gamma peak prominence (pooled)');

% Panel 4: burst prominence / coverage
nexttile(t,4);
yyaxis left
bar(T.burstProm);
ylabel('Median burst prominence');
yyaxis right
plot(T.burstCov,'o-','LineWidth',1.5);
ylabel('Burst coverage (QC breaths)');
set(gca,'XTick',1:nC,'XTickLabel',T.chanLabel,'XTickLabelRotation',45);
grid on; box off;
title('GammaBurst quality');

end



% ============================================================
% ============ FIGURES: GammaOsEvidence split ================
% ============================================================

function fig = plot_GammaOsEvidence_A1A2(chanDat, epochNames, gammaBandHz, baselineBandHz, excludeHzAroundPeak)

% -------------------- Inputs / defaults --------------------
if nargin < 5 || isempty(excludeHzAroundPeak), excludeHzAroundPeak = 5; end
if nargin < 4 || isempty(baselineBandHz), baselineBandHz = [25 58]; end
if nargin < 3 || isempty(gammaBandHz), gammaBandHz = [25 60]; end

nEpochs = numel(epochNames);

% QC breaths
if isfield(chanDat,'use') && ~isempty(chanDat.use)
    useVec = chanDat.use(:)==1;
else
    useVec = true(size(chanDat.fooof.gamma_peak_freq,1),1);
end
den = sum(useVec);

% Frequencies
if isfield(chanDat,'tf') && isfield(chanDat.tf,'frex')
    frex = double(chanDat.tf.frex(:));
else
    error('chanDat.tf.frex not found.');
end

% Peak frequencies per breath x epoch
G = [];

G = double(chanDat.fooof.gamma_peaks);
for ii = 1:size(G,1)
    for jj = 1:size(G,2)
        if isnan(G(ii,jj))
            G(ii,jj) = chanDat.fooof.gamma_peak_freq(ii,jj); 
        end
    end
end


% Flattened spectra in log10 space (breaths x epochs x frex)
flat_log10 = double(chanDat.fooof.spectra_flat_log10);


% -------------------- Compute spectral prominence per breath x epoch --------------------
% prom(b,e) in log10 units
[nBreaths, nEpochs2, nF] = size(flat_log10);
if nEpochs2 ~= nEpochs
    error('spectra_flat_log10 epoch dimension (%d) does not match gamma peaks (%d).', nEpochs2, nEpochs);
end
if numel(frex) ~= nF
    error('frex length (%d) does not match spectra_flat_log10 frex dim (%d).', numel(frex), nF);
end

prom = nan(nBreaths, nEpochs);     % spectral prominence (log10 units)
fiPk = nan(nBreaths, nEpochs);     % nearest frex index for the peak
pkHz = nan(nBreaths, nEpochs);     % snapped peak Hz (frex(fiPk))

baseBandMask = frex>=baselineBandHz(1) & frex<=baselineBandHz(2);
gammaMaskAll = frex>=gammaBandHz(1) & frex<=gammaBandHz(2);

for e = 1:nEpochs
    for b = 1:nBreaths
        if ~useVec(b), continue; end
        f0 = G(b,e);
        if ~isfinite(f0), continue; end

        % snap to frex
        [~,ii] = min(abs(frex - f0));
        fiPk(b,e) = ii;
        pkHz(b,e) = frex(ii);

        % baseline = median in baseline band excluding +-exclude around pk
        exclMask = frex >= (pkHz(b,e)-excludeHzAroundPeak) & frex <= (pkHz(b,e)+excludeHzAroundPeak);
        bmask = baseBandMask & ~exclMask;

        y = squeeze(flat_log10(b,e,:));
        if ~any(isfinite(y(gammaMaskAll))), continue; end

        base = median(y(bmask), 'omitnan');
        if ~isfinite(base), continue; end

        prom(b,e) = y(ii) - base;
    end
end

% -------------------- Figure layout (8 panels) --------------------
fig = figure('Color','w','Units','normalized','Position',[0.04 0.06 0.92 0.86]);
t = tiledlayout(fig, 2, 4, 'Padding','compact', 'TileSpacing','compact');

title(t, sprintf("GammaOsEvidence A1+A2 | %s | %s", string(chanDat.subID), chan_label(chanDat)), ...
    'Interpreter','none','FontWeight','bold');

% ---- (1) FOOOF prevalence bar (proportion + counts) ----
nexttile(t,1);
counts = zeros(1,nEpochs);
props  = nan(1,nEpochs);
for e=1:nEpochs
    counts(e) = sum(useVec & isfinite(chanDat.fooof.gamma_peaks(:,e)));
    props(e)  = counts(e)/max(den,1);
end
bar(props);
ylim([0 1]);
set(gca,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
ylabel('Proportion of breaths');
title('FOOOF gamma detection');
grid on; box off;
for e=1:nEpochs
    text(e, min(0.98, props(e)+0.04), sprintf('%d',counts(e)), ...
        'HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
end

% ---- (2) epoch x freq “histogram heatmap” FOOOF ----
nexttile(t,2);
fGmask = frex>=gammaBandHz(1) & frex<=gammaBandHz(2);
fG = frex(fGmask);
H = zeros(nEpochs, numel(fG));
for e=1:nEpochs
    x = chanDat.fooof.gamma_peaks(useVec & ...
        isfinite(chanDat.fooof.gamma_peaks(:,e)), e);
    for i=1:numel(x)
        [~,ii] = min(abs(fG - x(i)));
        H(e,ii) = H(e,ii) + 1;
    end
end
imagesc([], 1:nEpochs, H);

set(gca,'YTick',1:nEpochs,'YTickLabel',epochNames,'YDir','normal');
xticks(5:5:34)
xticklabels(round(fG(5:5:34)))
xlabel('Frequency (Hz)'); ylabel('Epoch');
title('Peak freq FOOOF');
if median(max(H, [], 2)) == 0
    e = 1; 
else
    e = median(max(H, [],2));
end
clim([0 e])
colorbar; box off;
hold on
maxFreqi = H*fG ./ sum(H,2);              % Hz mean
if sum(isnan(maxFreqi))==0
    maxFreqi = arrayfun(@(x) find(abs(fG - x)== ...
                                min(abs(fG - x)),1), maxFreqi);
    y = 1:nEpochs;
    plot(maxFreqi, y, 'k.-', 'LineWidth', 1);         % or 'k.-' if you want markers
end
hold off


% ---- (3) A1: epoch x freq “histogram heatmap” (QC peaks binned onto frex) ----
nexttile(t,3);
fGmask = frex>=gammaBandHz(1) & frex<=gammaBandHz(2);
fG = frex(fGmask);
H = zeros(nEpochs, numel(fG));
for e=1:nEpochs
     x = chanDat.fooof.gamma_peak_freq(useVec & ...
        isfinite(chanDat.fooof.gamma_peak_freq(:,e)), e);
    for i=1:numel(x)
        [~,ii] = min(abs(fG - x(i)));
        H(e,ii) = H(e,ii) + 1;
    end
end
imagesc([], 1:nEpochs, H);

set(gca,'YTick',1:nEpochs,'YTickLabel',epochNames,'YDir','normal');
xticks(5:5:34)
xticklabels(round(fG(5:5:34)))
xlabel('Frequency (Hz)'); ylabel('Epoch');
title('Peak freq MAX');
clim([0 min(max(H, [],2))])
colorbar; box off;
hold on
maxFreqi = H*fG ./ sum(H,2);              % Hz mean
if sum(isnan(maxFreqi))==0
    maxFreqi = arrayfun(@(x) find(abs(fG - x)== ...
                                min(abs(fG - x)),1), maxFreqi);
    y = 1:nEpochs;
    plot(maxFreqi, y, 'k.-', 'LineWidth', 1);         % or 'k.-' if you want markers
end
hold off



% ---- (4) epoch x freq broken down by fooof detected v. max detected peaks ----
nexttile(t,4);

% 
% counts = zeros(1,nEpochs);
% props  = nan(1,nEpochs);
% for e=1:nEpochs
%     counts(e) = sum(useVec & isfinite(chanDat.fooof.gamma_peaks(:,e)));
%     props(e)  = counts(e)/max(den,1);
% end

% ---- gather values per epoch ----
vals_F = cell(nEpochs,1);
vals_M = cell(nEpochs,1);
for e = 1:nEpochs
    curNanMask = isnan(chanDat.fooof.gamma_peaks(:,e)); 

    vals_F{e} = chanDat.fooof.gamma_peaks(...
                    useVec & isfinite(chanDat.fooof.gamma_peaks(:,e)),e);

    vals_M{e} = chanDat.fooof.gamma_peak_freq(useVec & curNanMask, e);
end

% colors: fooof detected v. max detected
colF = [0.20 0.45 0.80];
colM = [0.85 0.35 0.20];

% global x-range
allV = [];
for e=1:nEpochs
    allV = [allV; vals_F{e}(:); vals_M{e}(:)];
end
allV = allV(isfinite(allV));
xMin = prctile(allV, 5);
xMax = prctile(allV, 95);
if ~isfinite(xMin) || ~isfinite(xMax) || xMin==xMax
    xMin = min(allV); xMax = max(allV);
    if xMin==xMax, xMin = xMin-1; xMax = xMax+1; end
end

cla; hold on;

% positions: two boxes per epoch
xEpoch = 1:nEpochs;
dx = 0.18;      % horizontal offset within epoch
boxW = 0.28;    % box width
capW = 0.18;    % whisker cap width

% helper to draw one box at (x0) with color
for e = 1:nEpochs
    % --- Short ---
    v = vals_F{e};
    if numel(v) >= 5
        draw_box(v, xEpoch(e)-dx, colF, boxW, capW, xMin, xMax);
    end
    % --- Long ---
    v = vals_M{e};
    if numel(v) >= 5
        draw_box(v, xEpoch(e)+dx, colM, boxW, capW, xMin, xMax);
    end
end

    % axes formatting
    xlim([0.5 nEpochs+0.5]);
    ylim([xMin xMax]);
    set(gca,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
    ylabel('Frequency (Hz)');
    title('FOOOF v. Max peaks');
    grid on; box off;

    % legend (dummy patches)
    hS = patch(nan,nan,colF,'FaceAlpha',0.25,'EdgeColor',colF,'LineWidth',1.2);
    hL = patch(nan,nan,colM,'FaceAlpha',0.25,'EdgeColor',colM,'LineWidth',1.2);
    legend([hS hL], {'FOOOF','MAX'}, 'Location','best', 'Box','off');

    hold off;


% ---- (5) epoch x freq broken down by fooof detected v. max detected peaks ----
nexttile(t,5);

% 
% counts = zeros(1,nEpochs);
% props  = nan(1,nEpochs);
% for e=1:nEpochs
%     counts(e) = sum(useVec & isfinite(chanDat.fooof.gamma_peaks(:,e)));
%     props(e)  = counts(e)/max(den,1);
% end

% ---- gather values per epoch ----
vals_F = cell(nEpochs,1);
vals_M = cell(nEpochs,1);
for e = 1:nEpochs
    curNanMask = isnan(chanDat.fooof.gamma_peaks(:,e)); 

    vals_F{e} = prom(...
                    useVec & isfinite(chanDat.fooof.gamma_peaks(:,e)),e);

    vals_M{e} = prom(useVec & curNanMask, e);
end

% colors: fooof detected v. max detected
colF = [0.20 0.45 0.80];
colM = [0.85 0.35 0.20];

% global x-range
allV = [];
for e=1:nEpochs
    allV = [allV; vals_F{e}(:); vals_M{e}(:)];
end
allV = allV(isfinite(allV));
xMin = prctile(allV, 5);
xMax = prctile(allV, 95);
if ~isfinite(xMin) || ~isfinite(xMax) || xMin==xMax
    xMin = min(allV); xMax = max(allV);
    if xMin==xMax, xMin = xMin-1; xMax = xMax+1; end
end

cla; hold on;

% positions: two boxes per epoch
xEpoch = 1:nEpochs;
dx = 0.18;      % horizontal offset within epoch
boxW = 0.28;    % box width
capW = 0.18;    % whisker cap width

% helper to draw one box at (x0) with color
for e = 1:nEpochs
    % --- Short ---
    v = vals_F{e};
    if numel(v) >= 5
        draw_box(v, xEpoch(e)-dx, colF, boxW, capW, xMin, xMax);
    end
    % --- Long ---
    v = vals_M{e};
    if numel(v) >= 5
        draw_box(v, xEpoch(e)+dx, colM, boxW, capW, xMin, xMax);
    end
end

    % axes formatting
    xlim([0.5 nEpochs+0.5]);
    ylim([xMin xMax]);
    set(gca,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
    ylabel('Peak Prominence');
    title('FOOOF v. Max prominence');
    grid on; box off;

    % legend (dummy patches)
    hS = patch(nan,nan,colF,'FaceAlpha',0.25,'EdgeColor',colF,'LineWidth',1.2);
    hL = patch(nan,nan,colM,'FaceAlpha',0.25,'EdgeColor',colM,'LineWidth',1.2);
    legend([hS hL], {'FOOOF','MAX'}, 'Location','best', 'Box','off');

    hold off;


% ---- (6) FOOOF peak breaths: stacked heatmap ----
ax6 = nexttile(t,6); 

allSpect  = chanDat.fooof.spectra_flat_log10;   % breaths x epoch x freq
peakFooof = chanDat.fooof.gamma_peaks;          % breaths x epoch (FOOOF peaks)
peakMax   = G;                                  % breaths x epoch (MAX peaks)  <-- change if needed

% --- enforce desired epoch order if epochNames contain these labels (fallback = 1:nEpochs) ---
epochOrder = 1:nEpochs;
if exist('epochNames','var') && numel(epochNames)==nEpochs
    nm  = lower(string(epochNames));
    key = ["inhale rise","inhale fall","exhale rise","exhale fall","pause"];
    ord = [];
    for k = 1:numel(key)
        hit = find(contains(nm,key(k)), 1, 'first');
        if ~isempty(hit), ord(end+1) = hit; end 
    end
    if numel(ord)==nEpochs
        epochOrder = ord;
    end
end

% --- build stacked matrix for FOOOF ---
Mfoo = [];
boundsFoo = zeros(1,nEpochs);   % end-row index after each epoch block
labelsFoo = strings(1,nEpochs);
nRows = 0;

for ii = 1:nEpochs
    e = epochOrder(ii);
    labelsFoo(ii) = string(epochNames{e});

    keep = useVec & isfinite(peakFooof(:,e));             % only non-NaN breaths
    Se = squeeze(allSpect(keep, e, :));                   % breaths_kept x freq
    if isempty(Se)
        boundsFoo(ii) = nRows;
        continue
    end

    
  % Ensure both are [n x 300] (row = observation, col = 300 features)
    Mfoo = reshape(Mfoo, [], 300);
    Se   = reshape(Se,   [], 300);
    
    Mfoo = [Mfoo; Se];
    nRows = size(Mfoo, 1);
    boundsFoo(ii) = nRows;
end


imagesc([], 1:size(Mfoo,1), Mfoo);

xlabel('Frequency (Hz)');
ylabel('Breaths (stacked by epoch)');
title('FOOOF breath spectra'); 
box('off'); colorbar;

% dashed separators + epoch labels on y-axis (midpoints)
hold on;
for b = boundsFoo(1:end-1)
    if b > 0
        yline(b+0.5, 'k--', 'LineWidth', 1.5);
    end
end
hold off;

startsFoo = [1 boundsFoo(1:end-1)+1];
midsFoo   = (startsFoo + boundsFoo)/2;
validFoo  = boundsFoo > 0;
for b = 2:length(midsFoo)
    if midsFoo(b) == midsFoo(b-1)
        midsFoo(b) = midsFoo(b)+.001
    end
end
set(ax6,'YTick',midsFoo(validFoo),'YTickLabel',cellstr(labelsFoo(validFoo)));


% (optional) keep your original tick style
ix = 33:33:min(300,numel(frex));
xticks(ix);
xticklabels(round(frex(ix)));
xlim([125,300])

xRange6 = round(xlim(ax6));
colIdx6 = max(1,xRange6(1)) : min(size(Mfoo,2),xRange6(2));
valsFoo = Mfoo(:,colIdx6);
valsFoo = valsFoo(:);


% ---- (7) MAX peak breaths: stacked heatmap ----
ax7 = nexttile(t,7); 

allSpect  = chanDat.fooof.spectra_flat_log10;   % breaths x epoch x freq
peakFooof = chanDat.fooof.gamma_peaks;          % breaths x epoch (FOOOF peaks)
peakMax   = G;                                  % breaths x epoch (MAX peaks)  <-- change if needed

% --- enforce desired epoch order if epochNames contain these labels (fallback = 1:nEpochs) ---
epochOrder = 1:nEpochs;
if exist('epochNames','var') && numel(epochNames)==nEpochs
    nm  = lower(string(epochNames));
    key = ["inhale rise","inhale fall","exhale rise","exhale fall","pause"];
    ord = [];
    for k = 1:numel(key)
        hit = find(contains(nm,key(k)), 1, 'first');
        if ~isempty(hit), ord(end+1) = hit; end 
    end
    if numel(ord)==nEpochs
        epochOrder = ord;
    end
end

% --- build stacked matrix for MAX ---
Mfoo = [];
boundsFoo = zeros(1,nEpochs);   % end-row index after each epoch block
labelsFoo = strings(1,nEpochs);
nRows = 0;

for ii = 1:nEpochs
    e = epochOrder(ii);
    labelsFoo(ii) = string(epochNames{e});

    keep = useVec & isnan(peakFooof(:,e));             % only non-NaN breaths
    Se = squeeze(allSpect(keep, e, :));                   % breaths_kept x freq
    if isempty(Se)
        boundsFoo(ii) = nRows;
        continue
    end

    
   % Ensure both are [n x 300] (row = observation, col = 300 features)
    Mfoo = reshape(Mfoo, [], 300);
    Se   = reshape(Se,   [], 300);
    
    Mfoo = [Mfoo; Se];
    nRows = size(Mfoo, 1);
    boundsFoo(ii) = nRows;
end


imagesc([], 1:size(Mfoo,1), Mfoo);

xlabel('Frequency (Hz)');
ylabel('Breaths (stacked by epoch)');
title('MAX breath spectra'); 
box('off'); colorbar;

% dashed separators + epoch labels on y-axis (midpoints)
hold on;
for b = boundsFoo(1:end-1)
    if b > 0
        yline(b+0.5, 'k--', 'LineWidth', 1.5);
    end
end
hold off;

startsFoo = [1 boundsFoo(1:end-1)+1];
midsFoo   = (startsFoo + boundsFoo)/2;
validFoo  = boundsFoo > 0;
set(ax7, 'YTick',midsFoo(validFoo),'YTickLabel',cellstr(labelsFoo(validFoo)));

% (optional) keep your original tick style
ix = 33:33:min(300,numel(frex));
xticks(ix);
xticklabels(round(frex(ix)));
xlim([125,300])

xRange7 = round(xlim(ax7));
colIdx7 = max(1,xRange7(1)) : min(size(Mfoo,2),xRange7(2));
valsMax = Mfoo(:,colIdx7);
valsMax = valsMax(:);


allVals = [valsFoo(:); valsMax(:)];
cmin = prctile(allVals,5);
cmax = prctile(allVals,99);
if isfinite(cmin) && isfinite(cmax) && cmin < cmax
    set([ax6 ax7],'CLim',[cmin cmax]);
end

% ---- (8) Gamma timeseries: stacked heatmap ----
ax8 = nexttile(t,8);

% figure out common time-length (should be 50, but keep robust)
tN = min(size(chanDat.tf.powZ,2), size(chanDat.tf.breathSeg,2));

% two breath x time matrices (NaN init)
nBreaths = size(chanDat.tf.powZ,1);
MfooTS = nan(nBreaths, tN);
MmaxTS = nan(nBreaths, tN);

fooBreaths = [];
maxBreaths = [];

% gp  = chanDat.fooof.gamma_peaks;
% gpf = chanDat.fooof.gamma_peak_freq;

for b = 1:nBreaths
    
    % skip non-QC breaths
    if chanDat.use(b) ~= 1
        continue
    end
    
   peakHz = chanDat.gammaBurst.freqHz(b); 
   usedFooof = chanDat.gammaBurst.fooofBased(b); 
    % closest frequency index
    [~, fi] = min(abs(frex - peakHz));
    
    % gamma timecourse
    ts = squeeze(chanDat.tf.powZ(b, 1:tN, fi));
    
    % store into appropriate matrix
    if usedFooof
        MfooTS(b,:) = ts;
        fooBreaths(end+1) = b; 
    else
        MmaxTS(b,:) = ts;
        maxBreaths(end+1) = b; 
    end
end

% eliminate NaN rows
fooKeep = any(isfinite(MfooTS),2);
maxKeep = any(isfinite(MmaxTS),2);
Mfoo2 = MfooTS(fooKeep,:);
Mmax2 = MmaxTS(maxKeep,:);

% mean + SEM (timepoint-wise)
x = 1:tN;

fooN   = sum(isfinite(Mfoo2),1);
fooMu  = mean(Mfoo2,1,'omitnan');
fooSEM = std(Mfoo2,0,1,'omitnan') ./ sqrt(fooN);

maxN   = sum(isfinite(Mmax2),1);
maxMu  = mean(Mmax2,1,'omitnan');
maxSEM = std(Mmax2,0,1,'omitnan') ./ sqrt(maxN);

% plot mean gamma timeseries with SEM shading (left y-axis)
yyaxis(ax8,'left');
cla(ax8);
hold(ax8,'on');

fooIdx = isfinite(fooMu) & isfinite(fooSEM);
if any(fooIdx)
    xf = x(fooIdx);
    y1 = (fooMu - fooSEM); y2 = (fooMu + fooSEM);
    fill(ax8, [xf fliplr(xf)], [y1(fooIdx) fliplr(y2(fooIdx))], [0 0 0], ...
        'FaceAlpha',0.15,'EdgeColor','none');
    plot(ax8, xf, fooMu(fooIdx), '-', 'Color', [0 0 0], 'LineWidth', 2);
end

maxIdx = isfinite(maxMu) & isfinite(maxSEM);
if any(maxIdx)
    xm = x(maxIdx);
    y1 = (maxMu - maxSEM); y2 = (maxMu + maxSEM);
    fill(ax8, [xm fliplr(xm)], [y1(maxIdx) fliplr(y2(maxIdx))], [0.5 0.5 0.5], ...
        'FaceAlpha',0.20,'EdgeColor','none');
    plot(ax8, xm, maxMu(maxIdx), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
end

xlabel(ax8,'Time');
ylabel(ax8,'Gamma (z)');
title(ax8,'Gamma timeseries');
box(ax8,'off');

% mean respiration traces per set
if ~isempty(fooBreaths)
    fooBreaths = fooBreaths(:);
    meanRespFoo = mean(chanDat.tf.breathSeg(fooBreaths,1:tN), 1, 'omitnan');
else
    meanRespFoo = nan(1,tN);
end

if ~isempty(maxBreaths)
    maxBreaths = maxBreaths(:);
    meanRespMax = mean(chanDat.tf.breathSeg(maxBreaths,1:tN), 1, 'omitnan');
else
    meanRespMax = nan(1,tN);
end

% overlay respiration on right y-axis
yyaxis(ax8,'right');
plot(ax8, x, meanRespFoo, '--', 'Color', [0 0.6 0], 'LineWidth', 2);   % dashed green
plot(ax8, x, meanRespMax, ':',  'Color', [1 0.5 0], 'LineWidth', 2);   % dotted orange
ylabel(ax8,'Respiration');

hold(ax8,'off');


% ---- add legends for gamma + respiration ----
yyaxis(ax8,'left');  hold(ax8,'on');
hGammaFoo = plot(ax8, nan, nan, '-',  'Color', [0 0 0],       'LineWidth', 2);
hGammaMax = plot(ax8, nan, nan, '-',  'Color', [0.5 0.5 0.5], 'LineWidth', 2);

yyaxis(ax8,'right'); hold(ax8,'on');
hRespFoo  = plot(ax8, nan, nan, '--', 'Color', [0 0.6 0],     'LineWidth', 2);
hRespMax  = plot(ax8, nan, nan, ':',  'Color', [1 0.5 0],     'LineWidth', 2);

lg = legend(ax8, [hGammaFoo hGammaMax hRespFoo hRespMax], ...
    {'FOOOF (epochs 2/3)', 'Max', ...
     'Resp (FOOOF breaths)', 'Resp (MAX breaths)'}, ...
    'Location','best');
lg.Box = 'off';



% -------- local helper (nested or place as separate local function) --------
function draw_box(v, x0, col, boxW, capW, yMin, yMax)
    v = v(isfinite(v));
    if numel(v) < 5, return; end

    q1 = prctile(v,25);
    q2 = prctile(v,50);
    q3 = prctile(v,75);
    iqrV = q3-q1;

    loFence = q1 - 1.5*iqrV;
    hiFence = q3 + 1.5*iqrV;

    vC = v(v>=loFence & v<=hiFence);
    if isempty(vC)
        wLo = min(v); wHi = max(v);
    else
        wLo = min(vC); wHi = max(vC);
    end

    % clamp to plot range
    q1 = max(q1, yMin); q3 = min(q3, yMax);
    q2 = min(max(q2, yMin), yMax);
    wLo = max(wLo, yMin); wHi = min(wHi, yMax);

    % box
    xb = [x0-boxW/2 x0+boxW/2 x0+boxW/2 x0-boxW/2];
    yb = [q1 q1 q3 q3];
    patch(xb, yb, col, 'FaceAlpha',0.25, 'EdgeColor',col, 'LineWidth',1.2);

    % median
    plot([x0-boxW/2 x0+boxW/2], [q2 q2], 'Color',col, 'LineWidth',2.0);

    % whiskers
    plot([x0 x0], [wLo q1], 'Color',col, 'LineWidth',.5);
    plot([x0 x0], [q3 wHi], 'Color',col, 'LineWidth',.5);

    % caps
    plot([x0-capW/2 x0+capW/2], [wLo wLo], 'Color',col, 'LineWidth',.5);
    plot([x0-capW/2 x0+capW/2], [wHi wHi], 'Color',col, 'LineWidth',.5);

    % outliers (light jitter)
    out = v(v < loFence | v > hiFence);
    if ~isempty(out)
        jx = (rand(size(out))-0.5)*0.06;
        scatter(x0+jx, out, 8, col, 'filled', 'MarkerFaceAlpha',0.18, 'MarkerEdgeAlpha',0);
    end
end


% 
% 
% % ---- (4) A1: peak freq vs breath length scatter (colored by epoch) ----
% nexttile(t,4);
% if isfield(chanDat,'behDat') 
%     len = double(chanDat.behDat.length(:));
% else
%     len = nan(nBreaths,1);
% end
% hold on;
% cols = lines(nEpochs);
% for e=1:nEpochs
%     m = useVec & isfinite(G(:,e)) & isfinite(len);
%     scatter(len(m), G(m,e), 14, cols(e,:), 'filled', 'MarkerFaceAlpha',0.55);
% end
% hold off;
% xlabel('Breath length (s)');
% ylabel('Gamma peak freq (Hz)');
% legend(epochNames,'Location','bestoutside');
% title('A1. Peak freq vs breath length');
% grid on; box off;
% 
% % ---- (5) Prominence density by epoch (5 lines) ----
% nexttile(t,5);
% hold on;
% allProm = prom(useVec,:);
% if all(~isfinite(allProm(:)))
%     text(0.1,0.5,'No spectral prominence values available.','Units','normalized'); axis off;
% else
%     % grid for KDE
%     vAll = allProm(isfinite(allProm));
%     xMin = min(vAll); xMax = max(vAll);
%     if xMin==xMax, xMin = xMin-1; xMax = xMax+1; end
%     xGrid = linspace(xMin, xMax, 250);
% 
%     for e=1:nEpochs
%         v = prom(useVec & isfinite(prom(:,e)), e);
%         if numel(v) < 3, continue; end
%         [xGrid, y] = simple_kde(v, length(xGrid));
%         plot(xGrid, y, 'LineWidth', 1.8);
%     end
%     xlabel('Prominence (log10 units; peak minus baseline)');
%     ylabel('Density');
%     title('Prominence density by epoch');
%     legend(epochNames,'Location','southwest');
%     grid on; box off;
% end
% hold off;
% 
% % % ---- (2) A1: peak freq by epoch boxplot + MAD ----
% % nexttile(t,2);
% % X = []; grp = [];
% % for e=1:nEpochs
% %     x = G(useVec & isfinite(G(:,e)), e);
% %     X = [X; x(:)];
% %     grp = [grp; repmat(e, numel(x), 1)];
% % end
% % if isempty(X)
% %     text(0.1,0.5,'No gamma peaks found (QC breaths).','Units','normalized'); axis off;
% % else
% %     if exist('boxchart','file') == 2
% %         boxchart(grp, X);
% %         set(gca,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
% %     else
% %         boxplot(X, grp, 'Labels', cellstr(epochNames));
% %         set(gca,'XTickLabelRotation',25);
% %     end
% %     ylabel('Gamma peak frequency (Hz)');
% %     title('A1. Peak frequency by epoch');
% %     grid on; box off; hold on;
% %     for e=1:nEpochs
% %         x = G(useVec & isfinite(G(:,e)), e);
% %         md = robust_mad(x);
% %         yTop = max(x,[],'omitnan');
% %         % if isfinite(yTop)
% %         %     text(e, yTop + 1, sprintf('MAD=%.2f', md), ...
% %         %         'HorizontalAlignment','center','FontSize',9);
% %         % end
% %     end
% %     hold off;
% % end
% 
% 
% 
% 
% 
% % ---- (6) A2: Prominence boxplot by epoch ----
% nexttile(t,6);
% Xp = []; grpP = [];
% for e=1:nEpochs
%     v = prom(useVec & isfinite(prom(:,e)), e);
%     Xp = [Xp; v(:)];
%     grpP = [grpP; repmat(e, numel(v), 1)];
% end
% if isempty(Xp)
%     text(0.1,0.5,'No prominence values to boxplot.','Units','normalized'); axis off;
% else
%     if exist('boxchart','file') == 2
%         boxchart(grpP, Xp);
%         set(gca,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
%     else
%         boxplot(Xp, grpP, 'Labels', cellstr(epochNames));
%         set(gca,'XTickLabelRotation',25);
%     end
%     ylabel('Prominence (log10 units)');
%     title('A2. Prominence by epoch');
%     grid on; box off;
% end
% 

% 
% % ---- (8) A2: Primary vs secondary prominence comparison (time-domain burst prominence) ----
% nexttile(t,8);
% havePrim = isfield(chanDat,'gammaBurst') && isfield(chanDat.gammaBurst,'prominence');
% haveSec  = isfield(chanDat,'gammaBurstSecondary') && isfield(chanDat.gammaBurstSecondary,'prominence');
% 
% if ~(havePrim && haveSec)
%     text(0.08,0.55,'gammaBurst(.prominence) or gammaBurstSecondary(.prominence) missing.', 'Units','normalized');
%     axis off;
% else
%     p1 = double(chanDat.gammaBurst.prominence(:));
%     p2 = double(chanDat.gammaBurstSecondary.prominence(:));
%     m = useVec & isfinite(p1) & isfinite(p2);
% 
%     if ~any(m)
%         text(0.1,0.5,'No paired primary/secondary prominences (QC breaths).','Units','normalized'); axis off;
%     else
%         scatter(p1(m), p2(m), 16, 'filled', 'MarkerFaceAlpha',0.55);
%         hold on;
%         xl = xlim; yl = ylim;
%         lo = min([xl yl]); hi = max([xl yl]);
%         plot([lo hi],[lo hi],'k--','LineWidth',1.2);
%         xlim([lo hi]); ylim([lo hi]);
%         xlabel('Primary burst prominence (time-domain)');
%         ylabel('Secondary burst prominence (time-domain)');
%         title('A2. Primary vs secondary burst prominence');
%         grid on; box off;
%         hold off;
%     end
% end

end

function fig = plot_GammaOsEvidence_A3A4(chanDat, epochNames, gammaBandHz, baselineBandHz, excludeHzAroundPeak)

% -------------------- Inputs / defaults --------------------
if nargin < 5 || isempty(excludeHzAroundPeak), excludeHzAroundPeak = 5; end
if nargin < 4 || isempty(baselineBandHz), baselineBandHz = [25 58]; end
if nargin < 3 || isempty(gammaBandHz), gammaBandHz = [25 60]; end

nEpochs = numel(epochNames);

% QC breaths
if isfield(chanDat,'use') && ~isempty(chanDat.use)
    useVec = chanDat.use(:)==1;
else
    useVec = true(size(chanDat.fooof.gamma_peak_freq,1),1);
end
den = sum(useVec);

% Frequencies
if isfield(chanDat,'tf') && isfield(chanDat.tf,'frex')
    frex = double(chanDat.tf.frex(:));
else
    error('chanDat.tf.frex not found.');
end

% Breath lengths (for median split)

len = double(chanDat.behDat.length(:));


% -------------------- Peak freq matrices --------------------
% Detection matrix for *prevalence*:
% Prefer gamma_peaks presence (true "detected"); if missing entirely, fall back to gamma_peak_freq.
if isfield(chanDat,'fooof') && isfield(chanDat.fooof,'gamma_peaks') && ~isempty(chanDat.fooof.gamma_peaks)
    Gdet = double(chanDat.fooof.gamma_peaks);
elseif isfield(chanDat,'fooof') && isfield(chanDat.fooof,'gamma_peak_freq') && ~isempty(chanDat.fooof.gamma_peak_freq)
    Gdet = double(chanDat.fooof.gamma_peak_freq);
else
    error('Neither chanDat.fooof.gamma_peaks nor chanDat.fooof.gamma_peak_freq found.');
end

% Merged peak freq matrix for prominence / freq plots:
% elementwise prefer gamma_peaks, else gamma_peak_freq
G = double(Gdet); % start with whatever we used for detection
if isfield(chanDat,'fooof') && isfield(chanDat.fooof,'gamma_peaks') && ~isempty(chanDat.fooof.gamma_peaks) && ...
   isfield(chanDat,'fooof') && isfield(chanDat.fooof,'gamma_peak_freq') && ~isempty(chanDat.fooof.gamma_peak_freq)
    Gp = double(chanDat.fooof.gamma_peaks);
    Gf = double(chanDat.fooof.gamma_peak_freq);
    if isequal(size(Gp), size(Gf))
        m = ~isfinite(Gp);
        G = Gp;
        G(m) = Gf(m);
    else
        % sizes disagree: just keep Gdet and warn silently
        G = double(Gdet);
    end
end

% Flattened spectra in log10 space (breaths x epochs x frex)
if ~(isfield(chanDat,'fooof') && isfield(chanDat.fooof,'spectra_flat_log10') && ~isempty(chanDat.fooof.spectra_flat_log10))
    error('chanDat.fooof.spectra_flat_log10 not found.');
end
flat_log10 = double(chanDat.fooof.spectra_flat_log10);

% -------------------- Compute spectral prominence per breath x epoch --------------------
[nBreaths, nEpochs2, nF] = size(flat_log10);
if nEpochs2 ~= nEpochs
    error('spectra_flat_log10 epoch dimension (%d) does not match epochNames (%d).', nEpochs2, nEpochs);
end
if numel(frex) ~= nF
    error('frex length (%d) does not match spectra_flat_log10 frex dim (%d).', numel(frex), nF);
end

prom = nan(nBreaths, nEpochs);     % spectral prominence (log10 units)
fiPk = nan(nBreaths, nEpochs);     % nearest frex index for the peak
pkHz = nan(nBreaths, nEpochs);     % snapped peak Hz

baseBandMask = frex>=baselineBandHz(1) & frex<=baselineBandHz(2);
gammaMaskAll = frex>=gammaBandHz(1) & frex<=gammaBandHz(2);

for e = 1:nEpochs
    for b = 1:nBreaths
        if ~useVec(b), continue; end
        f0 = G(b,e);
        if ~isfinite(f0), continue; end

        [~,ii] = min(abs(frex - f0));
        fiPk(b,e) = ii;
        pkHz(b,e) = frex(ii);

        exclMask = frex >= (pkHz(b,e)-excludeHzAroundPeak) & frex <= (pkHz(b,e)+excludeHzAroundPeak);
        bmask = baseBandMask & ~exclMask;

        y = squeeze(flat_log10(b,e,:));
        if ~any(isfinite(y(gammaMaskAll))), continue; end

        base = median(y(bmask), 'omitnan');
        if ~isfinite(base), continue; end

        prom(b,e) = y(ii) - base;
    end
end

% -------------------- Condition vector + labels --------------------
[condVec, condLabels, uCond] = get_condition_vector(chanDat); 
hasCond = ~isempty(condVec) && any(isfinite(condVec));

% -------------------- Median split breath length (ignore condition) --------------------
mLen = median(len(useVec & isfinite(len)), 'omitnan');
isShort = useVec & isfinite(len) & (len <= mLen);
isLong  = useVec & isfinite(len) & (len >  mLen);

% -------------------- Figure layout (2 x 4) --------------------
fig = figure('Color','w','Units','normalized','Position',[0.04 0.06 0.94 0.86]);
t = tiledlayout(fig, 2, 4, 'Padding','compact', 'TileSpacing','compact');
title(t, sprintf("GammaOsEvidence A3+A4 | %s | %s", string(chanDat.subID), chan_label(chanDat)), ...
    'Interpreter','none','FontWeight','bold');

colsEpoch = lines(nEpochs);

% ============================================================
% Panel 1: Condition stratification prevalence (cond x epoch)
% ============================================================
nexttile(t,1);
if ~hasCond
    axis off;
    text(0,0.5,"Condition fields not found; skipping condition prevalence.",'FontSize',11);
else
    nCond = numel(uCond);
    P = nan(nCond, nEpochs);
    Ctot = zeros(nCond,1);
    Cpk  = zeros(nCond, nEpochs);

    for k = 1:nCond
        c = uCond(k);
        idxC = (condVec==c) & isfinite(condVec);
        Ctot(k) = sum(useVec & idxC);
        for e = 1:nEpochs
            Cpk(k,e) = sum(useVec & idxC & isfinite(Gdet(:,e)));
            P(k,e)   = Cpk(k,e) / max(Ctot(k),1);
        end
    end

    imagesc(1:nEpochs, 1:nCond, P);
    set(gca,'YDir','normal','XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
    set(gca,'YTick',1:nCond,'YTickLabel',cellstr(condLabels));
    xlabel('Epoch'); ylabel('Condition');
    title('A3. Prevalence by condition × epoch');
    caxis([0 1]); colorbar; box off;

    % annotate counts: "pk/total"
    for k = 1:nCond
        for e = 1:nEpochs
            txt = sprintf('%d/%d', Cpk(k,e), Ctot(k));
            text(e, k, txt, 'HorizontalAlignment','center', 'FontSize',8, 'Color','k');
        end
    end
end


% ============================================================
% Four panels: FOOOF vs non-FOOOF, epochs 2–3 vs outside 2–3
%   - Panels A/B: FOOOF-discovered peaks only
%   - Panels C/D: breaths missing FOOOF peaks in that epoch-set (fallback)
%   - Color-code points by condition (3 condensed conditions: matte orange/purple/green)
%   - Add per-condition best-fit lines
%   - Crop axes to middle 95% of points (per panel)
% ============================================================

E23   = [2 3];                     % inhale fall, exhale rise
Eout  = setdiff(1:nEpochs, E23);

% FOOOF availability mask (breaths x epochs)
fooAvail = isfinite(chanDat.fooof.gamma_peaks);

% Choose fallback peak/prom arrays (edit names here if yours differ)
if exist('pkHzMax','var') && exist('promMax','var')
    pkHzNF = pkHzMax;
    promNF = promMax;
elseif exist('pkHz_nonFooof','var') && exist('prom_nonFooof','var')
    pkHzNF = pkHz_nonFooof;
    promNF = prom_nonFooof;
else
    pkHzNF = pkHz;   % last-resort fallback
    promNF = prom;
end

% Condition colors (if present)
if exist('hasCond','var') && hasCond
    nCond = numel(condLabels);
    if ~exist('condCols','var') || size(condCols,1) ~= nCond
        condCols = lines(nCond);
    end
end

% Condense all conditions containing "shadow" into a single condition (in-place only)
condLabels = string(condLabels);
nCond = find(contains(lower(condLabels), "shadow"));

if ~isempty(nCond)
    condLabels(nCond(1)) = "shadow";
    condVec(ismember(condVec, nCond)) = nCond(1);

    if numel(nCond) > 1
        condVec = condVec - reshape(sum(condVec(:) > reshape(nCond(2:end), 1, []), 2), size(condVec));
        condLabels(nCond(2:end)) = [];
        condCols(nCond(2:end), :) = [];
    end
end

nCond = numel(condLabels);

condLabels(cellfun(@(x) contains(x, 'audio'), condLabels)) = 'audio'; 
condLabels(cellfun(@(x) contains(x, 'focus'), condLabels)) = 'focus'; 

% Set dark/matte colors for the (expected) three condensed conditions: orange, purple, green

    condCols =  [
        0.80 0.40 0.10  % dark matte orange
        0.45 0.20 0.60  % dark matte purple
        0.10 0.55 0.25  % dark matte green
    ];

% ----------------------------
% Panel A: FOOOF peaks in epochs 2–3
% ----------------------------
axA = nexttile(t);
hold(axA,'on');

xAll = []; yAll = [];
xByC = cell(nCond,1); yByC = cell(nCond,1);


for k = 1:nCond
    idxC = (condVec==k) & isfinite(condVec);
    x = []; 
    y = []; 
    for e = E23
        m = useVec & idxC & fooAvail(:,e) & isfinite(pkHz(:,e)) & isfinite(prom(:,e));
        x = [x; pkHz(m,e)]; 
        y = [y; prom(m,e)]; 
    end
    if length(x)>0
        scatter(axA, x, y, 28, condCols(k,:), 'filled', 'MarkerFaceAlpha',0.65);
        xAll = [xAll; x]; yAll = [yAll; y];
        xByC{k} = [xByC{k}; x]; yByC{k} = [yByC{k}; y];
    end
  
end


if ~isempty(xAll)
    xL = prctile(xAll,[2 98]); yL = prctile(yAll,[2 98]);
    if xL(1) == xL(2), xL = xL + [-1 1]*eps; end
    if yL(1) == yL(2), yL = yL + [-1 1]*eps; end
    xlim(axA, xL); ylim(axA, yL);
else
    xlim(axA,[gammaBandHz(1) gammaBandHz(2)]);
end

yL =  legend(axA, cellstr(condLabels), 'Location','bestoutside');
yL.AutoUpdate = 'off'; 

if exist('hasCond','var') && hasCond
    for k = 1:nCond
        if numel(xByC{k}) >= 2 && numel(unique(xByC{k})) >= 2
            p = polyfit(xByC{k}, yByC{k}, 1);
            xf = linspace(xL(1), xL(2), 100);
            plot(axA, xf, polyval(p,xf), '-', 'Color', condCols(k,:), 'LineWidth', 2);
        end
    end
end

xlabel(axA,'Peak frequency (Hz; snapped to frex)');
ylabel(axA,'Prominence (log10 units)');
title(axA,'FOOOF peaks (epochs 2–3)');
grid(axA,'on'); box(axA,'off');

hold(axA,'off');

% ----------------------------
% Panel B: FOOOF peaks outside epochs 2–3
% ----------------------------
axB = nexttile(t);
hold(axB,'on');

xAll = []; yAll = [];
xByC = cell(nCond,1); yByC = cell(nCond,1);

for k = 1:nCond
    idxC = (condVec==k) & isfinite(condVec);
    x = [];
    y = [];
    for e = Eout
        m = useVec & idxC & fooAvail(:,e) & isfinite(pkHz(:,e)) & isfinite(prom(:,e));
        x = [x; pkHz(m,e)];
        y = [y; prom(m,e)];
    end
    if length(x)>0
        scatter(axB, x, y, 28, condCols(k,:), 'filled', 'MarkerFaceAlpha',0.65);
        xAll = [xAll; x]; yAll = [yAll; y];
        xByC{k} = [xByC{k}; x]; yByC{k} = [yByC{k}; y];
    end
end

if ~isempty(xAll)
    xL = prctile(xAll,[2 98]); yL = prctile(yAll,[2 98]);
    if xL(1) == xL(2), xL = xL + [-1 1]*eps; end
    if yL(1) == yL(2), yL = yL + [-1 1]*eps; end
    xlim(axB, xL); ylim(axB, yL);
else
    xlim(axB,[gammaBandHz(1) gammaBandHz(2)]);
end

yL = legend(axB, cellstr(condLabels), 'Location','bestoutside');
yL.AutoUpdate = 'off';

if exist('hasCond','var') && hasCond
    for k = 1:nCond
        if numel(xByC{k}) >= 2 && numel(unique(xByC{k})) >= 2
            p = polyfit(xByC{k}, yByC{k}, 1);
            xf = linspace(xL(1), xL(2), 100);
            plot(axB, xf, polyval(p,xf), '-', 'Color', condCols(k,:), 'LineWidth', 2);
        end
    end
end

xlabel(axB,'Peak frequency (Hz; snapped to frex)');
ylabel(axB,'Prominence (log10 units)');
title(axB,'FOOOF peaks (outside epochs 2–3)');
grid(axB,'on'); box(axB,'off');

hold(axB,'off');

% ----------------------------
% Panel C: Non-FOOOF peaks for breaths with NO FOOOF peaks in epochs 2–3
% ----------------------------
axC = nexttile(t);
miss23 = useVec & ~any(fooAvail(:,E23), 2);
hold(axC,'on');

xAll = []; yAll = [];
xByC = cell(nCond,1); yByC = cell(nCond,1);

for k = 1:nCond
    idxC = (condVec==k) & isfinite(condVec);
    x = [];
    y = [];
    for e = E23
        m = miss23 & idxC & isfinite(pkHzNF(:,e)) & isfinite(promNF(:,e));
        x = [x; pkHzNF(m,e)];
        y = [y; promNF(m,e)];
    end
    if length(x)>0
        scatter(axC, x, y, 28, condCols(k,:), 'filled', 'MarkerFaceAlpha',0.65);
        xAll = [xAll; x]; yAll = [yAll; y];
        xByC{k} = [xByC{k}; x]; yByC{k} = [yByC{k}; y];
    end
end

if ~isempty(xAll)
    xL = prctile(xAll,[2 98]); yL = prctile(yAll,[2 98]);
    if xL(1) == xL(2), xL = xL + [-1 1]*eps; end
    if yL(1) == yL(2), yL = yL + [-1 1]*eps; end
    xlim(axC, xL); ylim(axC, yL);
else
    xlim(axC,[gammaBandHz(1) gammaBandHz(2)]);
end

yL = legend(axC, cellstr(condLabels), 'Location','bestoutside');
yL.AutoUpdate = 'off';

if exist('hasCond','var') && hasCond
    for k = 1:nCond
        if numel(xByC{k}) >= 2 && numel(unique(xByC{k})) >= 2
            p = polyfit(xByC{k}, yByC{k}, 1);
            xf = linspace(xL(1), xL(2), 100);
            plot(axC, xf, polyval(p,xf), '-', 'Color', condCols(k,:), 'LineWidth', 2);
        end
    end
end

xlabel(axC,'Peak frequency (Hz; snapped to frex)');
ylabel(axC,'Prominence (log10 units)');
title(axC,'Non-FOOOF peaks (breaths missing FOOOF in epochs 2–3)');
grid(axC,'on'); box(axC,'off');

hold(axC,'off');

% ----------------------------
% Panel D: Non-FOOOF peaks for breaths with NO FOOOF peaks outside epochs 2–3
% ----------------------------
axD = nexttile(t);
missOut = useVec & ~any(fooAvail(:,Eout), 2);
hold(axD,'on');

xAll = []; yAll = [];
xByC = cell(nCond,1); yByC = cell(nCond,1);

for k = 1:nCond
    idxC = (condVec==k) & isfinite(condVec);
    x = [];
    y = [];
    for e = Eout
        m = missOut & idxC & isfinite(pkHzNF(:,e)) & isfinite(promNF(:,e));
        x = [x; pkHzNF(m,e)];
        y = [y; promNF(m,e)];
    end
    if length(x)>0
        scatter(axD, x, y, 28, condCols(k,:), 'filled', 'MarkerFaceAlpha',0.65);
        xAll = [xAll; x]; yAll = [yAll; y];
        xByC{k} = [xByC{k}; x]; yByC{k} = [yByC{k}; y];
    end
end

if ~isempty(xAll)
    xL = prctile(xAll,[2 98]); yL = prctile(yAll,[2 98]);
    if xL(1) == xL(2), xL = xL + [-1 1]*eps; end
    if yL(1) == yL(2), yL = yL + [-1 1]*eps; end
    xlim(axD, xL); ylim(axD, yL);
else
    xlim(axD,[gammaBandHz(1) gammaBandHz(2)]);
end

yL = legend(axD, cellstr(condLabels), 'Location','bestoutside');
yL.AutoUpdate = 'off';

if exist('hasCond','var') && hasCond
    for k = 1:nCond
        if numel(xByC{k}) >= 2 && numel(unique(xByC{k})) >= 2
            p = polyfit(xByC{k}, yByC{k}, 1);
            xf = linspace(xL(1), xL(2), 100);
            plot(axD, xf, polyval(p,xf), '-', 'Color', condCols(k,:), 'LineWidth', 2);
        end
    end
end

xlabel(axD,'Peak frequency (Hz; snapped to frex)');
ylabel(axD,'Prominence (log10 units)');
title(axD,'Non-FOOOF peaks (breaths missing FOOOF outside epochs 2–3)');
grid(axD,'on'); box(axD,'off');

hold(axD,'off');





% ============================================================
% Panel 2: Ridge plot prominence by condition × epoch
% ============================================================
nexttile(t,6);
if ~hasCond
    axis off;
    text(0,0.5,"Condition fields not found; skipping condition × epoch prominence.",'FontSize',11);
else
    
    box_prom_epoch_condition(prom, useVec, condVec, uCond, condLabels, epochNames, colsEpoch);
    title('A3. Prominence by condition × epoch (horizontal boxplots)');
end


% ============================================================
% Panel 3: Median split prevalence by epoch (short vs long)
% ============================================================
nexttile(t,7);
propsS = nan(1,nEpochs); propsL = nan(1,nEpochs);
cntS = zeros(1,nEpochs); cntL = zeros(1,nEpochs);
denS = sum(isShort); denL = sum(isLong);

for e = 1:nEpochs
    cntS(e) = sum(isShort & isfinite(Gdet(:,e)));
    cntL(e) = sum(isLong  & isfinite(Gdet(:,e)));
    propsS(e) = cntS(e)/max(denS,1);
    propsL(e) = cntL(e)/max(denL,1);
end

B = bar([propsS(:) propsL(:)], 'grouped');
ylim([0 1]);
set(gca,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
ylabel('Proportion of breaths');
legend({'Short (<= median)','Long (> median)'}, 'Location','best');
title('A3. Prevalence: median-split breath length');
grid on; box off;

% annotate counts above bars
for e = 1:nEpochs
    x1 = B(1).XEndPoints(e); y1 = B(1).YEndPoints(e);
    x2 = B(2).XEndPoints(e); y2 = B(2).YEndPoints(e);
    text(x1, min(0.98,y1+0.04), sprintf('%d',cntS(e)), 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    text(x2, min(0.98,y2+0.04), sprintf('%d',cntL(e)), 'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end

% ============================================================
% Panel 4: Median split prominence by epoch (epoch-only labels; color = short/long)
% ============================================================
nexttile(t,8);

% ---- gather values per epoch ----
valsS = cell(nEpochs,1);
valsL = cell(nEpochs,1);
for e = 1:nEpochs
    valsS{e} = prom(isShort & isfinite(prom(:,e)), e);
    valsL{e} = prom(isLong  & isfinite(prom(:,e)), e);
end

% check any data
anyData = false;
for e = 1:nEpochs
    if numel(valsS{e})>=3 || numel(valsL{e})>=3
        anyData = true; break;
    end
end

if ~anyData
    axis off;
    text(0,0.5,"No prominence values for median split.",'FontSize',11);
else
    % colors: short/long
    colS = [0.20 0.45 0.80];
    colL = [0.85 0.35 0.20];

    % global x-range
    allV = [];
    for e=1:nEpochs
        allV = [allV; valsS{e}(:); valsL{e}(:)];
    end
    allV = allV(isfinite(allV));
    xMin = prctile(allV, 2);
    xMax = prctile(allV, 98);
    if ~isfinite(xMin) || ~isfinite(xMax) || xMin==xMax
        xMin = min(allV); xMax = max(allV);
        if xMin==xMax, xMin = xMin-1; xMax = xMax+1; end
    end

    cla; hold on;

    % positions: two boxes per epoch
    xEpoch = 1:nEpochs;
    dx = 0.18;      % horizontal offset within epoch
    boxW = 0.28;    % box width
    capW = 0.18;    % whisker cap width

    % helper to draw one box at (x0) with color
    for e = 1:nEpochs
        % --- Short ---
        v = valsS{e};
        if numel(v) >= 5
            draw_box(v, xEpoch(e)-dx, colS, boxW, capW, xMin, xMax);
        end
        % --- Long ---
        v = valsL{e};
        if numel(v) >= 5
            draw_box(v, xEpoch(e)+dx, colL, boxW, capW, xMin, xMax);
        end
    end

    % axes formatting
    xlim([0.5 nEpochs+0.5]);
    ylim([xMin xMax]);
    set(gca,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
    ylabel('Prominence (log10 units)');
    title('A3. Prominence: short vs long (color-coded; ignore condition)');
    grid on; box off;

    % legend (dummy patches)
    hS = patch(nan,nan,colS,'FaceAlpha',0.25,'EdgeColor',colS,'LineWidth',1.2);
    hL = patch(nan,nan,colL,'FaceAlpha',0.25,'EdgeColor',colL,'LineWidth',1.2);
    legend([hS hL], {'Short','Long'}, 'Location','best', 'Box','off');

    hold off;
end

% -------- local helper (nested or place as separate local function) --------
function draw_box(v, x0, col, boxW, capW, yMin, yMax)
    v = v(isfinite(v));
    if numel(v) < 5, return; end

    q1 = prctile(v,25);
    q2 = prctile(v,50);
    q3 = prctile(v,75);
    iqrV = q3-q1;

    loFence = q1 - 1.5*iqrV;
    hiFence = q3 + 1.5*iqrV;

    vC = v(v>=loFence & v<=hiFence);
    if isempty(vC)
        wLo = min(v); wHi = max(v);
    else
        wLo = min(vC); wHi = max(vC);
    end

    % clamp to plot range
    q1 = max(q1, yMin); q3 = min(q3, yMax);
    q2 = min(max(q2, yMin), yMax);
    wLo = max(wLo, yMin); wHi = min(wHi, yMax);

    % box
    xb = [x0-boxW/2 x0+boxW/2 x0+boxW/2 x0-boxW/2];
    yb = [q1 q1 q3 q3];
    patch(xb, yb, col, 'FaceAlpha',0.25, 'EdgeColor',col, 'LineWidth',1.2);

    % median
    plot([x0-boxW/2 x0+boxW/2], [q2 q2], 'Color',col, 'LineWidth',2.0);

    % whiskers
    plot([x0 x0], [wLo q1], 'Color',col, 'LineWidth',1.1);
    plot([x0 x0], [q3 wHi], 'Color',col, 'LineWidth',1.1);

    % caps
    plot([x0-capW/2 x0+capW/2], [wLo wLo], 'Color',col, 'LineWidth',1.1);
    plot([x0-capW/2 x0+capW/2], [wHi wHi], 'Color',col, 'LineWidth',1.1);

    % outliers (light jitter)
    out = v(v < loFence | v > hiFence);
    if ~isempty(out)
        jx = (rand(size(out))-0.5)*0.06;
        scatter(x0+jx, out, 8, col, 'filled', 'MarkerFaceAlpha',0.18, 'MarkerEdgeAlpha',0);
    end
end




% % ============================================================
% % Panel 5: A4 mean flattened spectra by epoch
% % ============================================================
% nexttile(t,5);
% mFlat = squeeze(mean(flat_log10(useVec,:,:), 1, 'omitnan')); % 5 x frex
% plot(frex, mFlat', 'LineWidth', 1.5);
% xlabel('Frequency (Hz)'); ylabel('Mean flattened log10(power)');
% legend(epochNames,'Location','best');
% title('A4. Mean flattened spectra (QC breaths)');
% grid on; box off;
% 
% % ============================================================
% % Panel 6: A4 aperiodic params by epoch (offset/exponent only; median ± MAD)
% %   - skips knee_freq
% %   - fixed y-lims: offset [4 8], exponent [1 2]
% %   - white background
% %   - larger panel + margins so labels/titles don't get cut off
% %   - bars colored by epoch (MATLAB default first 5 colors)
% % ============================================================
% axHost = nexttile(t,6);
% pos = axHost.Position;
% delete(axHost);
% 
% % White background panel
% pnl = uipanel(fig, 'Units','normalized', 'Position',pos, 'BorderType','none', 'BackgroundColor','w');
% 
% % Pull aperiodic params
% haveAP = isfield(chanDat,'fooof') && isfield(chanDat.fooof,'aperiodic_params') && ~isempty(chanDat.fooof.aperiodic_params);
% if ~haveAP
%     ax = axes('Parent',pnl,'Units','normalized','Position',[0 0 1 1], 'Color','w');
%     axis(ax,'off');
%     text(ax,0,0.5,'Aperiodic params missing: chanDat.fooof.aperiodic_params','FontSize',11,'Interpreter','none');
% else
%     ap = double(chanDat.fooof.aperiodic_params); % breaths x epochs x 3
%     off  = squeeze(ap(:,:,1));   % offset
%     expn = squeeze(ap(:,:,2));   % exponent
% 
%     % Robust summary: median ± MAD (QC breaths only)
%     stats = @(x) deal( ...
%         median(x(useVec & isfinite(x)),'omitnan'), ...
%         robust_mad(x(useVec & isfinite(x))) );
% 
%     % Epoch colors = MATLAB default order
%     colsEpoch = lines(max(nEpochs,5));
%     colsEpoch = colsEpoch(1:nEpochs,:);
% 
%     % Axes positions: give more margins so nothing clips
%     left = 0.16; width = 0.80;
%     ax1 = axes('Parent',pnl,'Units','normalized','Position',[left 0.55 width 0.40], 'Color','w');
%     ax2 = axes('Parent',pnl,'Units','normalized','Position',[left 0.10 width 0.40], 'Color','w');
% 
%     % ---------- OFFSET ----------
%     med = nan(1,nEpochs); md = nan(1,nEpochs);
%     for e=1:nEpochs, [med(e), md(e)] = stats(off(:,e)); end
% 
%     b1 = bar(ax1, med, 0.85); hold(ax1,'on');
%     % color each bar by epoch
%     for e=1:nEpochs
%         b1.CData(e,:) = colsEpoch(e,:);
%     end
%     errorbar(ax1, 1:nEpochs, med, md, '.k', 'LineWidth',1);
%     hold(ax1,'off');
% 
%     set(ax1,'XTick',1:nEpochs,'XTickLabel',[]);
%     ylabel(ax1,'offset');
%     title(ax1,'A4. Aperiodic params (median \pm MAD)','FontWeight','bold');
%     ylim(ax1,[4 8]);
%     grid(ax1,'on'); box(ax1,'off');
% 
%     % ---------- EXPONENT ----------
%     med(:)=nan; md(:)=nan;
%     for e=1:nEpochs, [med(e), md(e)] = stats(expn(:,e)); end
% 
%     b2 = bar(ax2, med, 0.85); hold(ax2,'on');
%     for e=1:nEpochs
%         b2.CData(e,:) = colsEpoch(e,:);
%     end
%     errorbar(ax2, 1:nEpochs, med, md, '.k', 'LineWidth',1);
%     hold(ax2,'off');
% 
%     set(ax2,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
%     ylabel(ax2,'exponent');
%     ylim(ax2,[1 2]);
%     grid(ax2,'on'); box(ax2,'off');
% end
% 
% 
% 
% % ============================================================
% % Panel 7: Peak freq distribution (short vs long) by epoch
% %   Style matches Panel 4: epoch labels only; color encodes short/long
% %   Assumes draw_box(v, x0, col, boxW, capW, yMin, yMax) is available
% % ============================================================
% nexttile(t,7);
% 
% % gather values per epoch
% valsS = cell(nEpochs,1);
% valsL = cell(nEpochs,1);
% for e = 1:nEpochs
%     valsS{e} = G(isShort & isfinite(G(:,e)), e);
%     valsL{e} = G(isLong  & isfinite(G(:,e)), e);
% end
% 
% % check any data
% anyData = false;
% for e = 1:nEpochs
%     if numel(valsS{e})>=3 || numel(valsL{e})>=3
%         anyData = true; break;
%     end
% end
% 
% if ~anyData
%     axis off;
%     text(0,0.5,"No peak-freq values for median split.",'FontSize',11);
% else
%     % colors: short/long (match Panel 4)
%     colS = [0.20 0.45 0.80];
%     colL = [0.85 0.35 0.20];
% 
%     % global y-range for readability
%     allV = [];
%     for e=1:nEpochs
%         allV = [allV; valsS{e}(:); valsL{e}(:)];
%     end
%     allV = allV(isfinite(allV));
%     if isempty(allV)
%         axis off;
%         text(0,0.5,"No finite peak-freq values.",'FontSize',11);
%         return
%     end
%     yMin = prctile(allV, 2);
%     yMax = prctile(allV, 98);
%     if ~isfinite(yMin) || ~isfinite(yMax) || yMin==yMax
%         yMin = min(allV); yMax = max(allV);
%         if yMin==yMax, yMin = yMin-1; yMax = yMax+1; end
%     end
% 
%     cla; hold on;
% 
%     % positions: two boxes per epoch
%     xEpoch = 1:nEpochs;
%     dx   = 0.18;
%     boxW = 0.28;
%     capW = 0.18;
% 
%     for e = 1:nEpochs
%         v = valsS{e};
%         if numel(v) >= 5
%             draw_box(v, xEpoch(e)-dx, colS, boxW, capW, yMin, yMax);
%         end
%         v = valsL{e};
%         if numel(v) >= 5
%             draw_box(v, xEpoch(e)+dx, colL, boxW, capW, yMin, yMax);
%         end
%     end
% 
%     xlim([0.5 nEpochs+0.5]);
%     ylim([yMin yMax]);
%     set(gca,'XTick',1:nEpochs,'XTickLabel',epochNames,'XTickLabelRotation',25);
%     ylabel('Peak frequency (Hz)');
%     title('A3. Peak freq: short vs long (color-coded)');
%     grid on; box off;
% 
%     % legend via dummy patches
%     hS = patch(nan,nan,colS,'FaceAlpha',0.25,'EdgeColor',colS,'LineWidth',1.2);
%     hL = patch(nan,nan,colL,'FaceAlpha',0.25,'EdgeColor',colL,'LineWidth',1.2);
%     legend([hS hL], {'Short','Long'}, 'Location','best', 'Box','off');
% 
%     hold off;
% end
% 
% 
% % ============================================================
% % Panel 8: Mean aperiodic fit (log10) by epoch (sanity)
% % ============================================================
% nexttile(t,8);
% if isfield(chanDat,'fooof') && isfield(chanDat.fooof,'aperiodic_fit') && ~isempty(chanDat.fooof.aperiodic_fit)
%     apFit = double(chanDat.fooof.aperiodic_fit); % breaths x epoch x frex (power)
%     apFitLog = log10(max(apFit, eps));
%     mAp = squeeze(mean(apFitLog(useVec,:,:), 1, 'omitnan')); % 5 x frex
%     plot(frex, mAp', 'LineWidth', 1.5);
%     xlabel('Frequency (Hz)'); ylabel('Mean log10(aperiodic fit)');
%     ylim([1 8])
%     legend(epochNames,'Location','best');
%     title('A4. Mean aperiodic fit (QC breaths)');
%     grid on; box off;
% else
%     axis off;
%     text(0,0.55,"chanDat.fooof.aperiodic_fit missing; can't plot aperiodic curves.",'FontSize',11);
% end
% 
% end
end
% -------------------- helpers --------------------

function box_prom_epoch_condition(prom, useVec, condVec, uCond, condLabels, epochNames, colsEpoch)
% BOX_PROM_EPOCH_CONDITION  Horizontal grouped boxplots: condition (rows) × epoch (colored boxes)
% BUT: merges all conditions whose label contains "shadow" into one aggregated "Shadow" condition.
%
% Inputs:
%   prom        [nBreaths x nEpochs]
%   useVec      [nBreaths x 1] QC mask
%   condVec     [nBreaths x 1] condition ID per breath
%   uCond       [nCond x 1] unique condition IDs (sorted)
%   condLabels  [nCond x 1] string/char labels for each uCond
%   epochNames  {1 x nEpochs} or string array
%   colsEpoch   [nEpochs x 3] colors (e.g., lines(nEpochs))

nEpochs = numel(epochNames);

% ---- normalize labels to string ----
condLabels = string(condLabels(:));
uCond = double(uCond(:));
condVec = double(condVec(:));
useVec = logical(useVec(:));

% ---- Build grouped conditions: merge any label containing "shadow" ----
isShadow = contains(lower(condLabels), "shadow");

shadowConds = uCond(isShadow);
nonShadowConds  = uCond(~isShadow);
nonShadowLabels = condLabels(~isShadow);

groups = cell(0,1);
groupLabels = string.empty(0,1);

% keep non-shadow as individual groups (in their original order)
for i = 1:numel(nonShadowConds)
    groups{end+1,1} = nonShadowConds(i); %#ok<AGROW>
    groupLabels(end+1,1) = nonShadowLabels(i); %#ok<AGROW>
end

% add one aggregated shadow group (if any)
if ~isempty(shadowConds)
    groups{end+1,1} = shadowConds;
    groupLabels(end+1,1) = "Shadow";
end

nGroups = numel(groups);
if nGroups == 0
    axis off;
    text(0,0.5,'No condition labels available.', 'FontSize',11);
    return
end

% --- collect global x-limits for readability ---
vAll = prom(useVec & isfinite(condVec), :);
vAll = vAll(isfinite(vAll));
if isempty(vAll)
    axis off;
    text(0,0.5,'No prominence values available.', 'FontSize',11);
    return
end
xMin = prctile(vAll, 2);
xMax = prctile(vAll, 98);
if ~isfinite(xMin) || ~isfinite(xMax) || xMin==xMax
    xMin = min(vAll); xMax = max(vAll);
    if xMin==xMax, xMin = xMin-1; xMax = xMax+1; end
end

hold on;

% Layout parameters within each group row
rowY = 1:nGroups;
dY = 0.07;                 % vertical offset step between epochs
boxH = 0.09;               % box half-height
whiskH = 0.05;             % whisker cap half-height

ctr = (nEpochs+1)/2;
epochOffsets = ((1:nEpochs) - ctr) * (2*dY);

% draw each group × epoch box
for k = 1:nGroups
    memberConds = groups{k};
    idxC = useVec & isfinite(condVec) & ismember(condVec, memberConds);

    for e = 1:nEpochs
        v = prom(idxC & isfinite(prom(:,e)), e);
        if numel(v) < 5
            continue
        end

        y0 = rowY(k) + epochOffsets(e);

        % box stats
        q1 = prctile(v,25);
        q2 = prctile(v,50);
        q3 = prctile(v,75);
        iqrV = q3 - q1;

        % Tukey whiskers
        loFence = q1 - 1.5*iqrV;
        hiFence = q3 + 1.5*iqrV;
        vClamped = v(v>=loFence & v<=hiFence);
        if isempty(vClamped)
            wLo = min(v);
            wHi = max(v);
        else
            wLo = min(vClamped);
            wHi = max(vClamped);
        end

        % clamp to plotting range
        q1 = max(q1, xMin); q3 = min(q3, xMax);
        q2 = min(max(q2, xMin), xMax);
        wLo = max(wLo, xMin); wHi = min(wHi, xMax);

        col = colsEpoch(e,:);

        % --- box (q1..q3) ---
        xBox = [q1 q3 q3 q1];
        yBox = [y0-boxH y0-boxH y0+boxH y0+boxH];
        patch(xBox, yBox, col, 'FaceAlpha',0.25, 'EdgeColor',col, 'LineWidth',1.0);

        % --- median line ---
        plot([q2 q2], [y0-boxH y0+boxH], 'Color',col, 'LineWidth',1.8);

        % --- whiskers ---
        plot([wLo q1], [y0 y0], 'Color',col, 'LineWidth',1.0);
        plot([q3 wHi], [y0 y0], 'Color',col, 'LineWidth',1.0);

        % --- whisker caps ---
        plot([wLo wLo], [y0-whiskH y0+whiskH], 'Color',col, 'LineWidth',1.0);
        plot([wHi wHi], [y0-whiskH y0+whiskH], 'Color',col, 'LineWidth',1.0);

        % optional outliers
        out = v(v < loFence | v > hiFence);
        out = out(isfinite(out));
        if ~isempty(out)
            jy = (rand(size(out)) - 0.5) * 0.06;
            scatter(out, y0+jy, 6, col, 'filled', 'MarkerFaceAlpha',0.25, 'MarkerEdgeAlpha',0);
        end
    end
end

% axis formatting
xlim([xMin xMax]);
ylim([0.5 nGroups+0.5]);
set(gca, 'YTick', rowY, 'YTickLabel', cellstr(groupLabels));
set(gca, 'TickDir','out');
xlabel('Prominence (log10 units; peak minus baseline)');
ylabel('Condition (shadow merged)');
grid on; box off;

% legend (epoch colors)
h = gobjects(nEpochs,1);
for e = 1:nEpochs
    h(e) = plot(nan,nan,'LineWidth',2,'Color',colsEpoch(e,:));
end
legend(h, epochNames, 'Location','best', 'Box','off');

hold off;

end


% ============================================================
% ============== FIGURES: BurstTiming split ==================
% ============================================================

function fig = plot_BurstTiming_B1B2(chanDat, epochNames)
useVec = chanDat.use(:)==1;
len = double(chanDat.behDat.length(:));

gb = chanDat.gammaBurst;
t0s = double(gb.t0_sec(:));
idxFull = double(gb.t0_idx_full(:));

% Analytic phase from lowRsp (FFT-based)
lowRsp = double(chanDat.trial.lowRsp); % breaths x time
phiAnal = sample_analytic_phase(lowRsp, idxFull);

% Template-index phase: nearest targIDX bin (1..50) -> angle
targIDX = double(chanDat.targIDX); % breaths x 50
phiTarg = targ_phase_from_targIDX(targIDX, idxFull);

% FOOOF present in epochs 2/3
gp = chanDat.fooof.gamma_peaks; % breaths x epochs
foo23 = chanDat.gammaBurst.fooofBased;

% QC breaths for timing plots
qc = useVec & isfinite(t0s) & isfinite(len);

% Condition vector + condensing (rely on your helper)
[condVec, condLabels, uCond] = get_condition_vector(chanDat); 
hasCond = ~isempty(condVec) && any(isfinite(condVec));

condLabels = string(condLabels);

% merge shadow/audio/focus into single labels each (and combine indices)
[condVec, condLabels] = merge_conditions_by_keyword(condVec, condLabels, "shadow");
[condVec, condLabels] = merge_conditions_by_keyword(condVec, condLabels, "audio");
[condVec, condLabels] = merge_conditions_by_keyword(condVec, condLabels, "focus");

nCond = numel(condLabels);

% Dark/matte colors for the (expected) three condensed conditions: orange, purple, green
condCols = [
    0.80 0.40 0.10  % dark matte orange
    0.45 0.20 0.60  % dark matte purple
    0.10 0.55 0.25  % dark matte green
];
if nCond > size(condCols,1)
    condCols = [condCols; lines(nCond-size(condCols,1))];
end

[nBreaths, nEpochs, nF] = size(chanDat.tf.spectra);

% prom = nan(nBreaths, nEpochs);     % spectral prominence (log10 units)
% fiPk = nan(nBreaths, nEpochs);     % nearest frex index for the peak
% pkHz = nan(nBreaths, nEpochs);     % snapped peak Hz
% frex = chanDat.tf.frex; 
% baseBandMask = frex>=baselineBandHz(1) & frex<=baselineBandHz(2);
% gammaMaskAll = frex>=gammaBandHz(1) & frex<=gammaBandHz(2);

% for e = 1:nEpochs
%     for b = 1:nBreaths
%         if ~useVec(b), continue; end
%         f0 = G(b,e);
%         if ~isfinite(f0), continue; end
% 
%         [~,ii] = min(abs(frex - f0));
%         fiPk(b,e) = ii;
%         pkHz(b,e) = frex(ii);
% 
%         exclMask = frex >= (pkHz(b,e)-excludeHzAroundPeak) & frex <= (pkHz(b,e)+excludeHzAroundPeak);
%         bmask = baseBandMask & ~exclMask;
% 
%         y = squeeze(flat_log10(b,e,:));
%         if ~any(isfinite(y(gammaMaskAll))), continue; end
% 
%         base = median(y(bmask), 'omitnan');
%         if ~isfinite(base), continue; end
% 
%         prom(b,e) = y(ii) - base;
%     end
% end







fig = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.9]);
t = tiledlayout(fig, 3, 3, 'Padding','compact', 'TileSpacing','compact');
title(t, sprintf("BurstTiming B1+B2 | %s | %s", chanDat.subID, chan_label(chanDat)), ...
    'Interpreter','none','FontWeight','bold');

% ============================================================
% Panel 1: Burst timing density by condition (shadow/audio/focus)
% ============================================================
ax1 = nexttile(t,1); hold(ax1,'on');
if any(qc)
    xGrid = linspace(prctile(t0s(qc),1), prctile(t0s(qc),99), 250);
else
    xGrid = linspace(0,1,250);
end

if hasCond
    for k = 1:nCond
        m = qc & isfinite(condVec) & (condVec==k);
        if nnz(m) >= 2
            f = ksdensity(t0s(m), xGrid);
            plot(ax1, xGrid, f, 'LineWidth', 2.5, 'Color', condCols(k,:));
        end
    end
    lg = legend(ax1, cellstr(condLabels), 'Location','bestoutside'); lg.AutoUpdate = 'off';
else
    if nnz(qc) >= 2
        f = ksdensity(t0s(qc), xGrid);
        plot(ax1, xGrid, f, 'LineWidth', 2.5, 'Color', [0 0 0]);
    end
end
xlabel(ax1,'Gamma burst time (s)'); ylabel(ax1,'Density');
title(ax1,'Burst timing density by condition (QC)');
grid(ax1,'on'); box(ax1,'off');
hold(ax1,'off');

% ============================================================
% Panel 2: Burst timing density by breath length (median split)
% ============================================================
ax2 = nexttile(t,2); hold(ax2,'on');
medLen = median(len(qc), 'omitnan');
isShort = qc & (len <= medLen);
isLong  = qc & (len >  medLen);

if nnz(isShort) >= 2
    f = ksdensity(t0s(isShort), xGrid);
    plot(ax2, xGrid, f, 'LineWidth', 2.5, 'Color', [0.15 0.15 0.15]);
end
if nnz(isLong) >= 2
    f = ksdensity(t0s(isLong), xGrid);
    plot(ax2, xGrid, f, 'LineWidth', 2.5, 'Color', [0.55 0.55 0.55]);
end
xlabel(ax2,'Gamma burst time (s)'); ylabel(ax2,'Density');
title(ax2,'Burst timing density: short vs long (median split)');
grid(ax2,'on'); box(ax2,'off');
lg = legend(ax2, {sprintf('Short (≤%.2fs)', medLen), sprintf('Long (>%.2fs)', medLen)}, 'Location','best'); 
lg.AutoUpdate = 'off';
hold(ax2,'off');

% ============================================================
% Panel 3: Burst timing density by FOOOF presence in epochs 2/3
% ============================================================
ax3 = nexttile(t,3); hold(ax3,'on');
hasFoo = qc & foo23;
noFoo  = qc & ~foo23;

if nnz(hasFoo) >= 2
    f = ksdensity(t0s(hasFoo), xGrid);
    plot(ax3, xGrid, f, 'LineWidth', 2.5, 'Color', [0 0.6 0]); % green
end
if nnz(noFoo) >= 2
    f = ksdensity(t0s(noFoo), xGrid);
    plot(ax3, xGrid, f, 'LineWidth', 2.5, 'Color', [0.75 0.45 0.10]); % orange-ish
end
xlabel(ax3,'Gamma burst time (s)'); ylabel(ax3,'Density');
title(ax3,'Burst timing density: FOOOF present vs absent');
grid(ax3,'on'); box(ax3,'off');
lg = legend(ax3, {sprintf('FOOOF present (n=%d)', nnz(hasFoo)), sprintf('FOOOF absent (n=%d)', nnz(noFoo))}, 'Location','best');
lg.AutoUpdate = 'off';
hold(ax3,'off');

% ============================================================
% Panel 4: Polar density (analytic phase) by condition
% ============================================================
ax = nexttile(t,4);
delete(ax);                 % remove the cartesian axes created by nexttile
p4 = polaraxes('Parent', t);
p4.Layout.Tile = 4;
hold(p4,'on');
if hasCond
    for k = 1:nCond
        m = qc & isfinite(condVec) & (condVec==k) & isfinite(phiAnal);
        if nnz(m) >= 5
            [th, rr] = polar_kde(phiAnal(m), 10);
            rr = rr ./ sum(rr); 
            polarplot(p4, th, rr, 'LineWidth', 2.0, 'Color', condCols(k,:));
        end
    end
    lg = legend(p4, cellstr(condLabels), 'Location','bestoutside'); lg.AutoUpdate = 'off';
else
    m = qc & isfinite(phiAnal);
    if nnz(m) >= 5
        [th, rr] = polar_kde(phiAnal(m), 256);
        polarplot(p4, th, rr, 'LineWidth', 2.0, 'Color', [0 0 0]);
    end
end
title(p4,'Analytic phase density @ burst');
hold(p4,'off');

% ============================================================
% Panel 5: Polar density (analytic phase) by FOOOF present/absent in epochs 2/3
% ============================================================
ax = nexttile(t,5);
delete(ax);                 % remove the cartesian axes created by nexttile
p5 = polaraxes('Parent', t);
p5.Layout.Tile = 5;
hold(p5,'on');

m1 = hasFoo & isfinite(phiAnal);
m2 = noFoo  & isfinite(phiAnal);

if nnz(m1) >= 5
    [th, rr] = polar_kde(phiAnal(m1), 10);
    rr = rr ./ sum(rr);
    polarplot(p5, th, rr, 'LineWidth', 2.0, 'Color', [0 0.6 0]); % present
end
if nnz(m2) >= 5
    [th, rr] = polar_kde(phiAnal(m2), 10);
    rr = rr ./ sum(rr);
    polarplot(p5, th, rr, 'LineWidth', 2.0, 'Color', [0.75 0.45 0.10]); % absent
end

title(p5,'Analytic phase density: FOOOF present vs absent');
lg = legend(p5, {sprintf('FOOOF present (n=%d)', nnz(m1)),...
    sprintf('FOOOF absent (n=%d)', nnz(m2))}, 'Location','westoutside');
lg.AutoUpdate = 'off';

hold(p5,'off');


% ============================================================
% Panel 6: Mean respiration waveform vs phase angle (deg): analytic vs breathSeg
% ============================================================
ax6 = nexttile(t,6); hold(ax6,'on');
degCenters = (0:49)/50 * 360;

% analytic-phase-sampled mean breath from lowRsp (per your indexing rule)
respAnalMat = nan(numel(useVec), 50);
for b = 1:numel(useVec)
    if useVec(b) ~= 1 || ~isfinite(len(b))
        continue
    end
    s0 = 1000;
    s1 = round(len(b)*500) + 1000;
    s1 = min(s1, size(lowRsp,2));
    if s1 <= s0
        continue
    end
    seg = lowRsp(b, s0:s1);
    if ~any(isfinite(seg))
        continue
    end
    ph = angle(hilbert(seg));
    deg = mod(ph * (180/pi), 360);
    respAnalMat(b,:) = mean_by_phase_bins(seg, deg, 50);
end
meanRespAnal = mean(respAnalMat(useVec,:), 1, 'omitnan');

% template-sampled mean breath from breathSeg (breaths x 50)
breathSeg = double(chanDat.tf.breathSeg(:,1:50));
meanRespSeg = mean(breathSeg(useVec,:), 1, 'omitnan');

plot(ax6, degCenters, meanRespAnal, '-', 'LineWidth', 2.5, 'Color', [0 0 0]);
plot(ax6, degCenters, meanRespSeg,  '-', 'LineWidth', 2.5, 'Color', [0.2 0.2 0.8]);
xlabel(ax6,'Phase angle (deg)');
ylabel(ax6,'Mean respiration (a.u.)');
title(ax6,'Mean respiration vs phase: analytic(lowRsp) vs breathSeg(50 bins)');
grid(ax6,'on'); box(ax6,'off');
xlim([0 360])
lg = legend(ax6, {'Analytic-phase sampled (lowRsp)', 'Template sampled (breathSeg)'}, 'Location','best');
lg.AutoUpdate = 'off';
hold(ax6,'off');

% ============================================================
% Panel 7: Polar density (targIDX phase) by condition
% ============================================================
ax = nexttile(t,7);
delete(ax);                 % remove the cartesian axes created by nexttile
p7 = polaraxes('Parent', t);
p7.Layout.Tile = 7;
hold(p7,'on');

if hasCond
    for k = 1:nCond
        m = qc & isfinite(condVec) & (condVec==k) & isfinite(phiTarg);
        if nnz(m) >= 5
            [th, rr] = polar_kde(phiTarg(m), 10);
            rr = rr ./ sum(rr);
            polarplot(p7, th, rr, 'LineWidth', 2.0, 'Color', condCols(k,:));
        end
    end
    lg = legend(p7, cellstr(condLabels), 'Location','bestoutside'); 
    lg.AutoUpdate = 'off';
else
    m = qc & isfinite(phiTarg);
    if nnz(m) >= 5
        [th, rr] = polar_kde(phiTarg(m), 10);
        rr = rr ./ sum(rr);
        polarplot(p7, th, rr, 'LineWidth', 2.0, 'Color', [0 0 0]);
    end
end
title(p7,'targIDX phase density @ burst');
hold(p7,'off');

% ============================================================
% Panel 8: Polar density (targIDX phase) by FOOOF present/absent in epochs 2/3
% ============================================================
ax = nexttile(t,8);
delete(ax);                 % remove the cartesian axes created by nexttile
p8 = polaraxes('Parent', t);
p8.Layout.Tile = 8;
hold(p8,'on');

m1 = hasFoo & isfinite(phiTarg);
m2 = noFoo  & isfinite(phiTarg);

if nnz(m1) >= 5
    [th, rr] = polar_kde(phiTarg(m1), 10);
    rr = rr ./ sum(rr);
    polarplot(p8, th, rr, 'LineWidth', 2.0, 'Color', [0 0.6 0]); % present
end
if nnz(m2) >= 5
    [th, rr] = polar_kde(phiTarg(m2), 10);
    rr = rr ./ sum(rr);
    polarplot(p8, th, rr, 'LineWidth', 2.0, 'Color', [0.75 0.45 0.10]); % absent
end

title(p8,'targIDX phase density: FOOOF present vs absent');
lg = legend(p8, {sprintf('FOOOF present (n=%d)', nnz(m1)), ...
    sprintf('FOOOF absent (n=%d)', nnz(m2))}, 'Location','westoutside');
lg.AutoUpdate = 'off';

hold(p8,'off');

% ============================================================
% Panel 9: Summary / notes
% ============================================================
ax9 = nexttile(t,9); axis(ax9,'off');
text(ax9, 0, 0.85, sprintf("QC breaths (timing): %d", nnz(qc)), 'FontSize',12);
text(ax9, 0, 0.70, sprintf("Median breath length: %.2fs", medLen), 'FontSize',12);
text(ax9, 0, 0.55, sprintf("FOOOF present in epochs 2/3: %d | absent: %d", nnz(hasFoo), nnz(noFoo)), 'FontSize',12);
text(ax9, 0, 0.35, "Panels 4/7: phase @ burst time; polar KDE lines", 'FontSize',10);
text(ax9, 0, 0.20, "Panel 6: mean respiration vs phase (0..360 deg), analytic(lowRsp) vs breathSeg", 'FontSize',10);
text(ax9, 0, 0.05, sprintf("Epochs: %s", strjoin(string(epochNames), ", ")), 'FontSize',9, 'Interpreter','none');

end

% ----------------- helpers (subfunctions) -----------------

function [condVec, condLabels] = merge_conditions_by_keyword(condVec, condLabels, kw)
condLabels = string(condLabels);
idx = find(contains(lower(condLabels), lower(string(kw))));
if isempty(idx)
    return
end
keep = idx(1);
condLabels(keep) = lower(string(kw));
condVec(ismember(condVec, idx)) = keep;
if numel(idx) > 1
    for d = idx(2:end)
        condVec(condVec > d) = condVec(condVec > d) - 1;
    end
    condLabels(idx(2:end)) = [];
end
end

function [theta, dens] = polar_kde(phi, nGrid)
%POLAR_KDE  (line-hist version)
%   [theta, dens] = polar_kde(phi, nGrid)
%   phi   : angles in radians (any range)
%   nGrid : number of bins
%   theta : bin-center angles (radians), length nGrid+1 (closed)
%   dens  : counts per bin,            length nGrid+1 (closed)

phi = phi(:);
phi = phi(isfinite(phi));

% define bin edges/centers on [0, 2*pi)
edges = linspace(0, 2*pi, nGrid+1);
centers = edges(1:end-1) + diff(edges)/2;

if isempty(phi)
    theta = [centers centers(1)];
    dens  = [nan(1,nGrid) nan];
    return
end

phi = mod(phi, 2*pi);

% histogram counts
dens = histcounts(phi, edges);

% close curve
theta = [centers centers(1)];
dens  = [dens dens(1)];
end


function y50 = mean_by_phase_bins(resp, deg, nBins)
% resp: 1xT, deg: 1xT in [0,360)
edges = linspace(0,360,nBins+1);
bin = discretize(deg, edges);
y50 = nan(1,nBins);
for j = 1:nBins
    m = (bin == j);
    if any(m)
        y50(j) = mean(resp(m), 'omitnan');
    end
end
end

% 
% 
% 
% function fig = plot_BurstTiming_B1B2(chanDat, epochNames)
% useVec = chanDat.use(:)==1;
% len = double(chanDat.behDat.length(:));
% 
% gb = chanDat.gammaBurst;
% t0s = double(gb.t0_sec(:));
% idxFull = double(gb.t0_idx_full(:));
% 
% % Analytic phase from lowRsp (FFT-based)
% lowRsp = double(chanDat.trial.lowRsp); % breaths x time
% phiAnal = sample_analytic_phase(lowRsp, idxFull);
% 
% % Template-index phase: nearest targIDX bin (1..50) -> angle
% targIDX = double(chanDat.targIDX); % breaths x 50
% phiTarg = targ_phase_from_targIDX(targIDX, idxFull);
% 
% % Median split by breath length (QC breaths)
% qc = useVec & isfinite(t0s);
% medLen = median(len(qc),'omitnan');
% isShort = qc & (len <= medLen);
% isLong  = qc & (len >  medLen);
% 
% fig = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.9]);
% t = tiledlayout(fig, 3, 3, 'Padding','compact', 'TileSpacing','compact');
% title(t, sprintf("BurstTiming B1+B2 | %s | %s", chanDat.subID, chan_label(chanDat)), 'Interpreter','none','FontWeight','bold');
% 
% % B1: timing histogram
% nexttile(t,1);
% histogram(t0s(qc), 30);
% xlabel('Gamma burst time (s)'); ylabel('Count');
% title('gammaBurst.t0\_sec (QC breaths)');
% grid on; box off;
% 
% % B1: t0 vs breath length scatter
% nexttile(t,2);
% scatter(len(qc), t0s(qc), 12, 'filled', 'MarkerFaceAlpha',0.6);
% xlabel('Breath length (s)'); ylabel('t0\_sec (s)');
% title('Burst timing vs breath length');
% grid on; box off;
% 
% % B1: timing proportion (optional)
% nexttile(t,3);
% prop = t0s ./ len;
% scatter(len(qc), prop(qc), 12, 'filled', 'MarkerFaceAlpha',0.6);
% xlabel('Breath length (s)'); ylabel('t0 / length');
% title('Burst timing as fraction of breath');
% grid on; box off;
% 
% % B2: polar analytic phase (short vs long)
% nexttile(t,4);
% polarhistogram(phiAnal(isShort), 18, 'Normalization','probability');
% title(sprintf("Analytic phase (short, n=%d)", sum(isShort)));
% 
% nexttile(t,5);
% polarhistogram(phiAnal(isLong), 18, 'Normalization','probability');
% title(sprintf("Analytic phase (long, n=%d)", sum(isLong)));
% 
% % B2: polar targIDX phase (short vs long)
% nexttile(t,7);
% polarhistogram(phiTarg(isShort), 18, 'Normalization','probability');
% title('targIDX-phase (short)');
% 
% nexttile(t,8);
% polarhistogram(phiTarg(isLong), 18, 'Normalization','probability');
% title('targIDX-phase (long)');
% 
% % B2: scatter comparing phase definitions (primary only, QC breaths)
% nexttile(t,6);
% scatter_wrap_phase(phiTarg(qc), phiAnal(qc));
% xlabel('targIDX-derived phase (rad)');
% ylabel('analytic phase (rad)');
% title('Phase definition comparison (QC breaths)');
% grid on; box off;
% 
% % Leave one tile for notes
% nexttile(t,9);
% axis off;
% text(0,0.8, sprintf("Median breath length=%.2fs", medLen), 'FontSize',12);
% text(0,0.6, "Analytic phase: angle(analytic(lowRsp))", 'FontSize',10);
% text(0,0.4, "targIDX phase: nearest of 50 phase indices -> angle", 'FontSize',10);
% 
% end
% 

function fig = plot_BurstTiming_B3B4(chanDat, epochNames, erpWinSec, bandHalfWidthHz)
useVec = chanDat.use(:)==1;

% Primary burst
gb = chanDat.gammaBurst;
t0Idx = double(gb.t0_idx_full(:));
f0    = double(gb.freqHz(:));
prom  = double(gb.prominence(:));

% Phase angles (use targIDX-derived here for robustness)
phiTarg = targ_phase_from_targIDX(double(chanDat.targIDX), t0Idx);
% phiTarg = targ_phase_from_targIDX(double(chanDat.targIDX), idxFull);
% B3: polar heatmap (phase vs prominence)
qc = useVec & isfinite(phiTarg) & isfinite(prom);
fig = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.9 0.9]);
t = tiledlayout(fig, 2, 2, 'Padding','compact', 'TileSpacing','compact');
title(t, sprintf("BurstTiming B3+B4 | %s | %s", chanDat.subID, chan_label(chanDat)), 'Interpreter','none','FontWeight','bold');

nexttile(t,1);


polar_heatmap(phiTarg(qc), prom(qc), 8, 5);
title('Polar density: phase × burst prominence (primary)');

% B4: ERP trough-locked (primary)
nexttile(t,2);
gbs = chanDat.gammaBurstSecondary;
t0Idx2 = double(gbs.t0_idx_full(:));
plot_dual_erp_raw(chanDat.trial.data, t0Idx, t0Idx2, useVec, erpWinSec, ...
    "Primary", "secondary");

% B4: ERP trough-locked (secondary)
nexttile(t,3)
fooofBased = gb.fooofBased; 
t0Idx = double(gb.t0_idx_full(:));
t0Idx2 = t0Idx; 
t0Idx(~fooofBased) = nan(sum(~fooofBased),1); 
t0Idx2(fooofBased) = NaN;
plot_dual_erp_raw(chanDat.trial.data, t0Idx, ...
    t0Idx2, useVec, erpWinSec, ...
    "fooof", "Max");

% B3: Primary vs secondary timing scatter (if available)
nexttile(t,4);
secondPhase = gb.phaseID==2;
t0Idx = double(gb.t0_idx_full(:));
t0Idx(~fooofBased) = nan; 
t0Idx2 = t0Idx; 
t0Idx(~secondPhase) = nan; 
t0Idx2(secondPhase) = nan;
plot_dual_erp_raw(chanDat.trial.data, t0Idx, ...
    t0Idx2, useVec, erpWinSec, ...
    "inhaleFall", "other");

end



% ============================================================
% ============= FIGURES: GammaEventStruct ====================
% ============================================================

function fig = plot_GammaEventStruct_Overall_LenSplit(chanDat, zPowCLim, itpcCLim)
useVec = chanDat.use(:)==1;
len = double(chanDat.behDat.length(:));
qc = useVec;

medLen = median(len(qc),'omitnan');
isShort = qc & (len <= medLen);
isLong  = qc & (len >  medLen);


[condVec, condLabels, uCond] = get_condition_vector(chanDat); %#ok<ASGLU>
hasCond = ~isempty(condVec) && any(isfinite(condVec));

condLabels = string(condLabels);

% merge shadow/audio/focus into single labels each (and combine indices)
[condVec, condLabels] = merge_conditions_by_keyword(condVec, condLabels, "shadow");
[condVec, condLabels] = merge_conditions_by_keyword(condVec, condLabels, "audio");
[condVec, condLabels] = merge_conditions_by_keyword(condVec, condLabels, "focus");

nCond = numel(condLabels);

% Dark/matte colors for the (expected) three condensed conditions: orange, purple, green
condCols = [
    0.80 0.40 0.10  % dark matte orange
    0.45 0.20 0.60  % dark matte purple
    0.10 0.55 0.25  % dark matte green
];


allNull = zeros(1000,length(qc),21, 100); 
parfor fi = 1:100
    fi
stem = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\';
baseName = [chanDat.sessID '_macro_breathingTask_' num2str(chanDat.chi)];
saveDirNull = fullfile(stem,'CHANDAT_processed','burstNullFiles');
nullFile = fullfile(saveDirNull, sprintf('%s_BurstNull_fi%03d.mat', baseName, chanDat.gammaLockTF.frexSelIdx(fi)));
test = load(nullFile);
test = test.nullOut.primary; 

allNull(:,:,:,fi) = test.phaseNull;
end

fig = figure('Color','w','Units','normalized','Position',[0.03 0.05 0.94 0.88]);
t = tiledlayout(fig, 4, 2, 'Padding','compact', 'TileSpacing','compact');
title(t, sprintf("GammaEventStruct | Overall + LenSplit | %s | %s", chanDat.subID, chan_label(chanDat)), 'Interpreter','none','FontWeight','bold');
qc = qc & chanDat.gammaBurst.fooofBased; 
frexMax = 25; 
frexSelect = chanDat.gammaLockTF.frexSel < frexMax; 
% ---------- Row 1: overall ----------
nexttile(t,1);
[zTF, itpcZ, fVec, tVec] = gammaLock_compute(chanDat, qc, allNull);
breathHz = 1 / mean(chanDat.behDat.length(qc));
[~, breathFi] = min(abs(chanDat.gammaLockTF.frexSel - breathHz));

imagesc(zTF(:, frexSelect)'); set(gca,'ydir','normal')
colorbar; clim([-2 2])
yticks(10:20:100); yticklabels(round(chanDat.gammaLockTF.frexSel(10:20:100)))
xticks(6:5:20);   xticklabels(-500:500:500)
yline(breathFi,'color','red','linestyle','--')
xlabel('time locked to gamma burst (ms)'); ylabel('Frequency (Hz)')
title('zPow primary (overall)');

nexttile(t,2);
imagesc(itpcZ(:, frexSelect)'); set(gca,'ydir','normal')
colorbar; clim([-2 2])
yticks(10:20:100); yticklabels(round(chanDat.gammaLockTF.frexSel(10:20:100)))
xticks(6:5:20);   xticklabels(-500:500:500)
yline(breathFi,'color','red','linestyle','--')
xlabel('time locked to gamma burst (ms)'); ylabel('Frequency (Hz)')
title('ITPC z (overall)');


% ---------- Row 2: cond 1 ----------
phase23 = condVec == 1;
keepB = qc & phase23;

% compute respiration line for this subset
breathHz = 1 / mean(chanDat.behDat.length(keepB));
[~, breathFi] = min(abs(chanDat.gammaLockTF.frexSel - breathHz));

nexttile(t,3);
[zTF, itpcZ, fVec, tVec] = gammaLock_compute(chanDat, keepB, allNull);

imagesc(zTF(:, frexSelect)'); set(gca,'ydir','normal')
colorbar; clim([-2 2])
yticks(10:20:100); yticklabels(round(chanDat.gammaLockTF.frexSel(10:20:100)))
xticks(6:5:20);   xticklabels(-500:500:500)
yline(breathFi,'color','red','linestyle','--')
xlabel('time locked to gamma burst (ms)'); ylabel('Frequency (Hz)')
title('zPow primary (cond 1)');

nexttile(t,4);
imagesc(itpcZ(:, frexSelect)'); set(gca,'ydir','normal')
colorbar;  clim([-2 2])
yticks(10:20:100); yticklabels(round(chanDat.gammaLockTF.frexSel(10:20:100)))
xticks(6:5:20);   xticklabels(-500:500:500)
yline(breathFi,'color','red','linestyle','--')
xlabel('time locked to gamma burst (ms)'); ylabel('Frequency (Hz)')
title('ITPC z (cond 1)');


% ---------- Row 3: cond 2 ----------
phase23 = condVec == 2;
keepB = qc & phase23;

breathHz = 1 / mean(chanDat.behDat.length(keepB));
[~, breathFi] = min(abs(chanDat.gammaLockTF.frexSel - breathHz));

nexttile(t,5);
[zTF, itpcZ, fVec, tVec] = gammaLock_compute(chanDat, keepB, allNull);

imagesc(zTF(:, frexSelect)'); set(gca,'ydir','normal')
colorbar; clim([-2 2])
yticks(10:20:100); yticklabels(round(chanDat.gammaLockTF.frexSel(10:20:100)))
xticks(6:5:20);   xticklabels(-500:500:500)
yline(breathFi,'color','red','linestyle','--')
xlabel('time locked to gamma burst (ms)'); ylabel('Frequency (Hz)')
title('zPow primary (cond 2)');

nexttile(t,6);
imagesc(itpcZ(:, frexSelect)'); set(gca,'ydir','normal')
colorbar;  clim([-2 2])
yticks(10:20:100); yticklabels(round(chanDat.gammaLockTF.frexSel(10:20:100)))
xticks(6:5:20);   xticklabels(-500:500:500)
yline(breathFi,'color','red','linestyle','--')
xlabel('time locked to gamma burst (ms)'); ylabel('Frequency (Hz)')
title('ITPC z (cond 2)');


% ---------- Row 4: cond 3 ----------
phase23 = condVec == 3;
keepB = qc & phase23;

breathHz = 1 / mean(chanDat.behDat.length(keepB));
[~, breathFi] = min(abs(chanDat.gammaLockTF.frexSel - breathHz));

nexttile(t,7);
[zTF, itpcZ, fVec, tVec] = gammaLock_compute(chanDat, keepB, allNull);

imagesc(zTF(:, frexSelect)'); set(gca,'ydir','normal')
colorbar; clim([-2 2])
yticks(10:20:100); yticklabels(round(chanDat.gammaLockTF.frexSel(10:20:100)))
xticks(6:5:20);   xticklabels(-500:500:500)
yline(breathFi,'color','red','linestyle','--')
xlabel('time locked to gamma burst (ms)'); ylabel('Frequency (Hz)')
title('zPow primary (cond 3)');

nexttile(t,8);
imagesc(itpcZ(:, frexSelect)'); set(gca,'ydir','normal')
colorbar; clim([-2 2])
yticks(10:20:100); yticklabels(round(chanDat.gammaLockTF.frexSel(10:20:100)))
xticks(6:5:20);   xticklabels(-500:500:500)
yline(breathFi,'color','red','linestyle','--')
xlabel('time locked to gamma burst (ms)'); ylabel('Frequency (Hz)')
title('ITPC z (cond 3)');






end


function fig = plot_GammaEventStruct_QualitySplit(chanDat, zPowCLim, itpcCLim)
useVec = chanDat.use(:)==1;

gb = chanDat.gammaBurst;
snr = nan(numel(useVec),1);
if isfield(gb,'snr')
    snr = double(gb.snr(:));
else
    snr = double(gb.prominence(:)); % fallback
end
qc = useVec & isfinite(snr);
medQ = median(snr(qc),'omitnan');
isHi = qc & (snr >= medQ);
isLo = qc & (snr <  medQ);

[zHi, itpcHi, fVec, tVec] = gammaLock_compute(chanDat, isHi);
[zLo, itpcLo] = gammaLock_compute(chanDat, isLo);

fig = figure('Color','w','Units','normalized','Position',[0.03 0.05 0.94 0.88]);
t = tiledlayout(fig, 2, 3, 'Padding','compact', 'TileSpacing','compact');
title(t, sprintf("GammaEventStruct | QualitySplit | %s | %s", chanDat.subID, chan_label(chanDat)), 'Interpreter','none','FontWeight','bold');

nexttile(t,1); imagesc(tVec, fVec, zHi'); set(gca,'YDir','normal');
caxis(zPowCLim); colorbar; xlabel('Time (s)'); ylabel('Freq (Hz)'); title(sprintf('zPow (highQ, n=%d)', sum(isHi)));
nexttile(t,2); imagesc(tVec, fVec, itpcHi'); set(gca,'YDir','normal');
caxis(itpcCLim); colorbar; xlabel('Time (s)'); ylabel('Freq (Hz)'); title('ITPC diff (highQ)');
nexttile(t,3);
[~,i0] = min(abs(tVec-0));
m = zHi(i0,:); s = stderr_over_breaths(chanDat, isHi, i0);
hold on; plot(fVec,m,'LineWidth',1.5); fill_between(fVec,m-s,m+s,0.2); hold off;
xlabel('Freq (Hz)'); ylabel('Mean zPow ± SEM'); title('zPow profile at t=0 (highQ)');
grid on; box off;

nexttile(t,4); imagesc(tVec, fVec, zLo'); set(gca,'YDir','normal');
caxis(zPowCLim); colorbar; xlabel('Time (s)'); ylabel('Freq (Hz)'); title(sprintf('zPow (lowQ, n=%d)', sum(isLo)));
nexttile(t,5); imagesc(tVec, fVec, itpcLo'); set(gca,'YDir','normal');
caxis(itpcCLim); colorbar; xlabel('Time (s)'); ylabel('Freq (Hz)'); title('ITPC diff (lowQ)');
nexttile(t,6);
m = zLo(i0,:); s = stderr_over_breaths(chanDat, isLo, i0);
hold on; plot(fVec,m,'LineWidth',1.5); fill_between(fVec,m-s,m+s,0.2); hold off;
xlabel('Freq (Hz)'); ylabel('Mean zPow ± SEM'); title('zPow profile at t=0 (lowQ)');
grid on; box off;

end



% ============================================================
% ====================== COMPUTATIONS ========================
% ============================================================

function [promEpoch, promPooled] = compute_gamma_prominence(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
% Returns median prominence per epoch and pooled median across all epochs
[promByEpoch, promAll] = compute_gamma_prominence_full(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak);
nEpochs = numel(promByEpoch);

promEpoch = nan(1,nEpochs);
for e=1:nEpochs
    promEpoch(e) = median(promByEpoch{e}, 'omitnan');
end
promPooled = median(promAll, 'omitnan');
end


function [promByEpoch, promAll, promRaw] = compute_gamma_prominence_full(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
% Returns per-epoch prominence samples (cell), pooled sample vector, and raw per-breath×epoch matrix

Gdet = double(chanDat.fooof.gamma_peaks);        % NaN if no peak
Gmax = double(chanDat.fooof.gamma_peak_freq);    % fallback
flat = double(chanDat.fooof.spectra_flat_log10); % breaths x epochs x frex
frex = double(chanDat.tf.frex(:));

nBreaths = size(Gdet,1);
nEpochs  = size(Gdet,2);

% Choose peak frequency per breath/epoch: detected if possible else max
Guse = Gdet;
missing = ~isfinite(Guse);
Guse(missing) = Gmax(missing);

% Precompute masks for baseline band
maskBaseBand = frex>=baselineBandHz(1) & frex<=baselineBandHz(2);

promRaw = nan(nBreaths, nEpochs);
for b=1:nBreaths
    if ~useVec(b), continue; end
    for e=1:nEpochs
        fpk = Guse(b,e);
        if ~isfinite(fpk), continue; end
        [~, fiPk] = min(abs(frex - fpk));

        % baseline mask excludes +/- excludeHzAroundPeak around peak
        maskEx = maskBaseBand & ~(frex >= (fpk-excludeHzAroundPeak) & frex <= (fpk+excludeHzAroundPeak));

        base = median(flat(b,e,maskEx), 3, 'omitnan');
        if ~isfinite(base), continue; end

        promRaw(b,e) = flat(b,e,fiPk) - base;
    end
end

promByEpoch = cell(1,nEpochs);
promAll = [];
for e=1:nEpochs
    x = promRaw(useVec, e);
    x = x(isfinite(x));
    promByEpoch{e} = x(:);
    promAll = [promAll; x(:)];
end

end


function [zTF, itpcZ, fVec, tVec] = gammaLock_compute(chanDat, keepBreaths, allNull)
%GAMMALOCK_COMPUTE
%   zTF   : primary power z-scored vs secondary (mean across kept breaths) [time x freq]
%   itpcZ : ITPC of PRIMARY phases z-scored vs null ITPC distribution from allNull [time x freq]
%
% Inputs
%   chanDat     : struct with chanDat.gammaLockTF fields
%   keepBreaths : logical vector (nBreaths x 1) selecting breaths to use
%   allNull     : [nShuf x nBreaths x time x freq] phase values (radians),
%                 where nShuf=1000 permutation shuffles
%
% Notes
%   - allNull is assumed to already be phase values for the *null* alignment.
%   - ITPC null is computed across breaths for each shuffle, yielding 1000 ITPC maps.

GL = chanDat.gammaLockTF;

powP = double(GL.pow_primary);
powS = double(GL.pow_secondary);
phP  = double(GL.phase_primary);

keepBreaths = keepBreaths(:);
idx = find(keepBreaths);

% --- basic checks ---
nBreaths = size(powP,1);
if numel(keepBreaths) ~= nBreaths
    error('keepBreaths must have length size(GL.pow_primary,1).');
end

if ndims(allNull) ~= 4
    error('allNull must be 4-D: [nShuf x nBreaths x time x freq].');
end
if size(allNull,2) ~= nBreaths
    error('allNull second dimension (breaths) must match size(GL.pow_primary,1).');
end
if size(allNull,3) ~= size(phP,2) || size(allNull,4) ~= size(phP,3)
    error('allNull [time x freq] dims must match GL.phase_primary.');
end

% Restrict to subset of breaths
powP = powP(idx,:,:);
powS = powS(idx,:,:);
phP  = phP(idx,:,:);

allNull = double(allNull(:, idx, :, :));  % [nShuf x nKept x time x freq]

% -------------------------
% Power zTF (unchanged logic)
% -------------------------
muS = squeeze(mean(powS, 1, 'omitnan'));    % [time x freq]
sdS = squeeze(std(powS, 0, 1, 'omitnan'));  % [time x freq]
sdS(sdS==0 | ~isfinite(sdS)) = NaN;

% z for primary breaths, then average across breaths
zBreath = (powP - permute(muS,[3 1 2])) ./ permute(sdS,[3 1 2]); % [breath x time x freq]
zTF     = squeeze(mean(zBreath, 1, 'omitnan'));                  % [time x freq]

% -------------------------
% Observed ITPC (PRIMARY)
% -------------------------
Zp       = exp(1i * phP);                           % [breath x time x freq]
itpcObs  = abs(squeeze(mean(Zp, 1, 'omitnan')));    % [time x freq]

% -------------------------
% Null ITPC distribution from allNull
% -------------------------
Znull     = exp(1i * allNull);                      % [nShuf x breath x time x freq]
itpcNull  = abs(squeeze(mean(Znull, 2, 'omitnan')));% [nShuf x time x freq]

muNull = squeeze(mean(itpcNull, 1, 'omitnan'));     % [time x freq]
sdNull = squeeze(std( itpcNull, 0, 1, 'omitnan'));  % [time x freq]
sdNull(sdNull==0 | ~isfinite(sdNull)) = NaN;

% Z-score ITPC relative to null
itpcZ = (itpcObs - muNull) ./ sdNull;               % [time x freq]

% Axes
fVec = double(GL.frexSel(:));
tVec = double(GL.offsets_out(:)) ./ double(GL.fs_out); % seconds

end



function s = stderr_over_breaths(chanDat, keepBreaths, timeIndex)
% SEM across breaths at a given time index for zPow (computed like gammaLock_compute)
GL = chanDat.gammaLockTF;
powP = double(GL.pow_primary);
powS = double(GL.pow_secondary);

idx = find(keepBreaths(:));
powP = powP(idx,:,:);
powS = powS(idx,:,:);

muS = squeeze(mean(powS, 1, 'omitnan'));  % [time x freq]
sdS = squeeze(std(powS, 0, 1, 'omitnan'));% [time x freq]
sdS(sdS==0 | ~isfinite(sdS)) = NaN;

% z at chosen time index
mu = muS(timeIndex,:);
sd = sdS(timeIndex,:);
z = squeeze((powP(:,timeIndex,:) - mu) ./ sd); % [breath x freq]

% SEM
n = sum(isfinite(z),1);
s = std(z,0,1,'omitnan') ./ sqrt(max(n,1));
end



% ============================================================
% ====================== PLOTTING UTILS ======================
% ============================================================

function plot_dual_erp_raw(trialData, idxAFull, idxBFull, useVec, erpWinSec, labelA, labelB)
%PLOT_DUAL_ERP_RAW  Plot two raw ERPs overlaid, aligned to provided indices.
%
% Inputs
%   trialData : breaths x time (raw)
%   idxAFull  : nBreaths x 1 indices for condition A alignment (no trough search)
%   idxBFull  : nBreaths x 1 indices for condition B alignment (no trough search)
%   useVec    : nBreaths x 1 logical; breaths to include (e.g., chanDat.use==1)
%   erpWinSec : half-window in seconds (ERP spans [-erpWinSec, +erpWinSec])
%   labelA    : string for legend (optional)
%   labelB    : string for legend (optional)
%
% Notes
%   - Does NOT do trough alignment.
%   - Uses idxAFull / idxBFull exactly as provided.
%   - Assumes your fixed trial window is -2..10 s (12 s total) and infers fs = nT/12.

if nargin < 6 || isempty(labelA), labelA = 'A'; end
if nargin < 7 || isempty(labelB), labelB = 'B'; end

X = double(trialData);                 % breaths x time
nBreaths = size(X,1);
nT       = size(X,2);

idxAFull = idxAFull(:);
idxBFull = idxBFull(:);
useVec   = useVec(:);

if numel(idxAFull) ~= nBreaths || numel(idxBFull) ~= nBreaths || numel(useVec) ~= nBreaths
    error('idxAFull, idxBFull, and useVec must all have length nBreaths = size(trialData,1).');
end

% Infer fs from fixed -2..10 s window (12 s total)
fs = round(nT/12);

winSamp = round(erpWinSec * fs);
tVec    = (-winSamp:winSamp) / fs;

% Collect epochs (raw only)
segA = nan(nBreaths, numel(tVec));
segB = nan(nBreaths, numel(tVec));

for b = 1:nBreaths
    if ~useVec(b), continue; end

    % --- Condition A ---
    idx0 = idxAFull(b);
    if isfinite(idx0) && idx0 >= 1 && idx0 <= nT
        lo = idx0 - winSamp;
        hi = idx0 + winSamp;
        if lo >= 1 && hi <= nT
            segA(b,:) = X(b, lo:hi);
        end
    end

    % --- Condition B ---
    idx0 = idxBFull(b);
    if isfinite(idx0) && idx0 >= 1 && idx0 <= nT
        lo = idx0 - winSamp;
        hi = idx0 + winSamp;
        if lo >= 1 && hi <= nT
            segB(b,:) = X(b, lo:hi);
        end
    end
end

% Mean ± SEM (computed per timepoint using available breaths)
mA = mean(segA, 1, 'omitnan');
mB = mean(segB, 1, 'omitnan');

nA = sum(isfinite(segA), 1);
nB = sum(isfinite(segB), 1);

sA = std(segA, 0, 1, 'omitnan') ./ sqrt(max(nA,1));
sB = std(segB, 0, 1, 'omitnan') ./ sqrt(max(nB,1));

% Plot
ax = gca;
cla(ax);
hold(ax,'on');

p1 = plot(tVec, mA, 'LineWidth', 1.5);
fill_between(tVec, mA - sA, mA + sA, 0.20);

p2 = plot(tVec, mB, 'LineWidth', 1.5);
fill_between(tVec, mB - sB, mB + sB, 0.20);

xlabel('Time from provided index (s)');
ylabel(['Raw (uV)']);
title(sprintf('Raw ERP overlay: %s vs %s', labelA, labelB));
grid on; box off;
legend([p1 p2], {sprintf('%s (n~%d)', labelA, max(nA)), sprintf('%s (n~%d)', labelB, max(nB))}, ...
    'Location','best');
hold(ax,'off');

end



function plot_trough_locked_erp(trialData, t0IdxFull, f0, useVec, erpWinSec, bandHalfWidthHz, labelStr)
% Align each breath to nearest trough of bandpassed signal within ±half-cycle of burst index.
% Then plot mean ± SEM for bandpassed and raw signals (dual y-axis).

X = double(trialData); % breaths x time
nBreaths = size(X,1);
nT = size(X,2);

% Estimate fs from trial time axis if available? Not passed. Assume 500? NO:
% Infer fs from typical length if tim exists elsewhere is safer; here we infer by 12s window:
% Better: use chanDat.tim but not passed. We'll estimate by nT/12 seconds (since -2..10).
fs = round(nT/12); % robust given your fixed window

halfCycleSamp = round((fs ./ max(f0,eps)) * 0.5);

winSamp = round(erpWinSec * fs);
tVec = (-winSamp:winSamp) / fs;

% Collect epochs
rawSeg = nan(nBreaths, numel(tVec));
bpSeg  = nan(nBreaths, numel(tVec));

for b=1:nBreaths
    if ~useVec(b), continue; end
    idx0 = t0IdxFull(b);
    if ~isfinite(idx0) || idx0 < 1 || idx0 > nT, continue; end
    fb = f0(b);
    if ~isfinite(fb) || fb<=0, continue; end

    % bandpass around fb ± bandHalfWidthHz (FFT-based)
    x = X(b,:);
    xbp = bandpass_fft(x, fs, fb-bandHalfWidthHz, fb+bandHalfWidthHz);

    hc = halfCycleSamp(b);
    lo = max(1, idx0 - hc);
    hi = min(nT, idx0 + hc);

    [~, rel] = min(xbp(lo:hi)); % nearest trough
    idxT = lo + rel - 1;

    segLo = idxT - winSamp;
    segHi = idxT + winSamp;
    if segLo < 1 || segHi > nT, continue; end

    rawSeg(b,:) = x(segLo:segHi);
    bpSeg(b,:)  = xbp(segLo:segHi);
end

% Mean ± SEM
mRaw = mean(rawSeg, 1,'omitnan');
sRaw = std(rawSeg, 0, 1,'omitnan') ./ sqrt(max(sum(isfinite(rawSeg),1),1));

mBP  = mean(bpSeg, 1,'omitnan');
sBP  = std(bpSeg, 0, 1,'omitnan') ./ sqrt(max(sum(isfinite(bpSeg),1),1));

yyaxis left
hold on;
plot(tVec, mBP, 'LineWidth',1.5);
fill_between(tVec, mBP-sBP, mBP+sBP, 0.2);
ylabel('Bandpassed (a.u.)');
yyaxis right
plot(tVec, mRaw, 'LineWidth',1.0);
fill_between(tVec, mRaw-sRaw, mRaw+sRaw, 0.15);
ylabel('Raw (a.u.)');
hold off;

xlabel('Time from nearest trough (s)');
title(sprintf('Trough-locked ERP (%s)', labelStr));
grid on; box off;

end


function polar_heatmap(theta, r, nThetaBins, nRBins)
%POLAR_HEATMAP  Polar "heatmap" (2D density) in (theta,r).
%   theta: angles (rad)
%   r    : radius values (e.g., prominence)
%   nThetaBins, nRBins: number of bins

if nargin < 3 || isempty(nThetaBins), nThetaBins = 36; end
if nargin < 4 || isempty(nRBins),     nRBins     = 20; end

theta = theta(:);
r     = r(:);

ok = isfinite(theta) & isfinite(r);
theta = theta(ok);
r     = r(ok);

if isempty(theta)
    axis off; text(0,0.5,"No data.",'FontSize',12); return;
end

% Wrap to [-pi, pi)
theta = mod(theta + pi, 2*pi) - pi;

% --- bin edges ---
tEdges = linspace(-pi, pi, nThetaBins+1);

% (Your snippet had rHi commented out; this keeps intent sane.)
rLo = min(r);
rHi = prctile(r, 90);
if rLo == rHi, rHi = rLo + 1; end

rMin = min(0, rLo);
rMax = rHi;
rEdges = linspace(rMin, rMax, nRBins+1);

% Keep your "nudge off the edges" behavior
r(r <= rLo) = rLo + 1e-5;
r(r >= rHi) = rHi - 1e-5;

% --- YOUR H calculation (unchanged) ---
H = zeros(nThetaBins, nRBins);
for tc = 1:nThetaBins
    for rc = 1:nRBins
        H(tc, rc) = sum(theta>=tEdges(tc) & theta<tEdges(tc+1) ...
            & r>=rEdges(rc) & r<rEdges(rc+1));
    end
end
H = H.'; % (nRBins x nThetaBins)

% =========================
% FIX: plot using EDGES so we get nRBins radial faces (not nRBins-1)
% =========================

% Build vertex grid from edges: size = (nRBins+1) x (nThetaBins+1)
[TT, RR] = meshgrid(tEdges, rEdges);
XX = RR .* cos(TT);
YY = RR .* sin(TT);

% Map bin counts to face colors by storing H at the lower-left vertex of each face.
% Create CData same size as XX/YY; last row/col are just padding (not used as lower-left of any face).
C = zeros(size(XX));
C(1:end-1, 1:end-1) = H;

% Pad last row/col (values won't affect face colors in flat mode, but keeps C fully defined)
C(end, 1:end-1) = H(end, :);
C(1:end-1, end) = H(:, end);
C(end, end)     = H(end, end);

ax = gca;
cla(ax);
hold(ax,'on');

surf(ax, XX, YY, zeros(size(XX)), C, 'EdgeColor','none', 'FaceColor','flat');
view(ax, 2);
axis(ax, 'equal');
ax.XLim = [-rMax, rMax];
ax.YLim = [-rMax, rMax];
ax.Box  = 'off';
ax.XTick = [];
ax.YTick = [];

colormap(ax, parula);
cb = colorbar(ax);
cb.Label.String = 'Counts';
title(ax, 'Density (counts)');

% --- FIX #2: add polar grid + angle/radius labels --- 
% Angle ticks (degrees) and labels 
angDeg = 0:45:315; 
angRad = deg2rad(angDeg); 
labRad = 1.08 * rMax; 
for k = 1:numel(angRad) 
    th = angRad(k); 
    % radial spokes 
    plot(ax, [0, rMax*cos(th)], [0, rMax*sin(th)], 'k:', ... 
        'HandleVisibility','off'); 
    % angle labels 
    text(ax, labRad*cos(th), labRad*sin(th), sprintf('%d°', angDeg(k)), ... 
        'HorizontalAlignment','center', 'VerticalAlignment','middle', ... 
        'FontSize', 10, 'Color', 'k', 'Clipping','off'); 
end


hold(ax,'off');
end




function fill_between(x, y1, y2, alphaVal)
x = x(:); y1 = y1(:); y2 = y2(:);
patch([x; flipud(x)], [y1; flipud(y2)], 'k', 'FaceAlpha',alphaVal, 'EdgeColor','none');
end


function scatter_wrap_phase(x, y)
% Wrap y into [-pi,pi) and show scatter with reference lines
x = wrapToPi(x(:));
y = wrapToPi(y(:));
scatter(x,y,12,'filled','MarkerFaceAlpha',0.5);
xlim([-pi pi]); ylim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]); yticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
line([-pi pi],[0 0],'Color',[0.7 0.7 0.7]);
line([0 0],[-pi pi],'Color',[0.7 0.7 0.7]);
end



% ============================================================
% ======================== SIGNAL UTILS ======================
% ============================================================

function y = bandpass_fft(x, fs, fLo, fHi)
% Simple FFT bandpass for real vector x
x = x(:);
n = numel(x);
X = fft(x);

f = (0:n-1)'*(fs/n);
mask = (f>=fLo & f<=fHi) | (f>=fs-fHi & f<=fs-fLo);
X(~mask) = 0;

y = real(ifft(X));
y = y(:).';
end


function phi = sample_analytic_phase(xMat, idxVec)
% xMat: breaths x time
% idxVec: breaths x 1 indices (1-based) into columns
% returns phase at idx for each breath

[nB,nT] = size(xMat);
phi = nan(nB,1);
idxVec = round(idxVec);

for b=1:nB
    idx = idxVec(b);
    if ~isfinite(idx) || idx<1 || idx>nT, continue; end
    x = xMat(b,:);
    % analytic signal via FFT Hilbert method
    z = analytic_signal_fft(x);
    phi(b) = angle(z(idx));
end

end


function z = analytic_signal_fft(x)
% Analytic signal (Hilbert) via FFT, no toolbox
x = x(:);
n = numel(x);
X = fft(x);

H = zeros(n,1);
if mod(n,2)==0
    H(1) = 1; H(n/2+1) = 1;
    H(2:n/2) = 2;
else
    H(1) = 1;
    H(2:(n+1)/2) = 2;
end

z = ifft(X .* H);
end


function phi = targ_phase_from_targIDX(targIDX, idxFull)
% targIDX: breaths x 50 (raw indices)
% idxFull: breaths x 1 (raw index)
% returns angle in radians from nearest targIDX column index

nB = size(targIDX,1);
phi = nan(nB,1);
for b=1:nB
    idx = idxFull(b);
    if ~isfinite(idx), continue; end
    v = targIDX(b,:);
    v = v(:);
    ok = isfinite(v);
    if ~any(ok), continue; end
    [~,k] = min(abs(v(ok) - idx));
    % map nearest among ok indices back to original 1..50
    okIdx = find(ok);
    k50 = okIdx(k); % 1..50
    % convert 1..50 -> angle
    phi(b) = (double(k50)-1)/50 * 2*pi;
end
phi = wrapToPi(phi);
end



% ============================================================
% ==================== CONDITION UTILITIES ===================
% ============================================================

function col = getBehCol(behDat, varName, nBreaths)
col = nan(nBreaths,1);
try
    if istable(behDat) && any(strcmp(behDat.Properties.VariableNames,varName))
        tmp = behDat.(varName);
        col = double(tmp(:));
    elseif isstruct(behDat) && isfield(behDat,varName)
        tmp = behDat.(varName);
        col = double(tmp(:));
    end
catch
end
col = padOrTruncate(col, nBreaths, NaN);
end

function A = padOrTruncate(A, n, fillVal)
A = A(:);
if numel(A) < n
    A = [A; fillVal*ones(n-numel(A),1)];
elseif numel(A) > n
    A = A(1:n);
end
end

function s = safeStr(x)
try
    if isstring(x), s = char(x); return; end
    if ischar(x), s = x; return; end
    if iscell(x) && ~isempty(x), s = safeStr(x{1}); return; end
    if isnumeric(x) && isscalar(x), s = num2str(x); return; end
    s = char(string(x));
catch
    s = 'NA';
end
if isempty(s), s = 'NA'; end
end

function col = getBehColCell(behDat, varName, nBreaths)
col = cell(nBreaths,1); [col{:}] = deal('NA');
try
    if istable(behDat) && any(strcmp(behDat.Properties.VariableNames,varName))
        tmp = behDat.(varName);
    elseif isstruct(behDat) && isfield(behDat,varName)
        tmp = behDat.(varName);
    else
        tmp = [];
    end
    if isempty(tmp), return; end

    tmp = tmp(:);
    for i = 1:min(nBreaths,numel(tmp))
        if iscell(tmp)
            col{i} = safeStr(tmp{i});
        else
            col{i} = safeStr(tmp(i));
        end
    end
catch
end
end

function [condVec, condLabels, uCond] = get_condition_vector(chanDat)

    if ~strcmp(chanDat.task, 'breathingTask')
        nBreaths = size(chanDat.trial.data, 1);
        task = getBehColCell(chanDat.behDat,'sniffLabel',nBreaths);
        cond = getBehCol(chanDat.behDat,'condition',nBreaths);
        if sum(isnan(cond))==length(cond)
            cond = getBehColCell(chanDat.behDat, 'sniffLabel', nBreaths); 
            [~,~,cond] = unique(cond, 'stable'); 
        end
        uCond = unique(cond(isfinite(cond)));
        uCond = uCond(:)';
        nCond = numel(uCond);
        condVec = nan(nBreaths,1); 
        useVec = chanDat.use==1; 

        countsPass = zeros(1,nCond);
        condLabels = strings(1,nCond);

        for k = 1:nCond
            c = uCond(k);
            idx = (cond==c) & isfinite(cond);
            condVec(idx) = c; 
            countsPass(k) = sum(useVec & idx);
            j = find(idx,1,'first');
            
            curTask = task{j}; 
            
            if ~isempty(j)
                condLabels(k) = [num2str(c) ...
                    curTask];
            else
                condLabels(k) = sprintf('c%d', round(c));
            end
        end

    else
    nBreaths = size(chanDat.trial.data, 1);
    warp = getBehColCell(chanDat.behDat,'warp',nBreaths);
    cond = getBehCol(chanDat.behDat,'condition',nBreaths);
    task = getBehColCell(chanDat.behDat,'task',nBreaths);
    nmth = getBehColCell(chanDat.behDat,'noseMouth',nBreaths);
    uCond = unique(cond(isfinite(cond)));
    uCond = uCond(:)';
    nCond = numel(uCond);
    condVec = nan(nBreaths,1); 
    useVec = chanDat.use==1; 

    countsPass = zeros(1,nCond);
    condLabels = strings(1,nCond);

        for k = 1:nCond
            c = uCond(k);
            idx = (cond==c) & isfinite(cond);
            condVec(idx) = c; 
            countsPass(k) = sum(useVec & idx);
            j = find(idx,1,'first');

            if isnan(warp{j})
                curWarp = ''; 
            else
                curWarp = warp{j}; 
            end
            if strcmp(task{j}, 'shadow')
                if strcmp(chanDat.behDat.shadowFile{j}, 'focusedResp')
                    curTask = 'FocShadow';
                else
                    curTask = 'AudShadow';
                end
            else
                curTask = task{j}; 
            end
            if ~isempty(j)
                condLabels(k) = [num2str(c) ...
                    curTask nmth{j}, ...
                   curWarp];
            else
                condLabels(k) = sprintf('c%d', round(c));
            end
        end

    end
end

% ---------- helpers ----------
function c = to_cellstr(v, n)
% Convert string/char/cellstr to cellstr length n (best effort)
if isstring(v)
    c = cellstr(v(:));
elseif ischar(v)
    c = cellstr(string(v));
elseif iscell(v)
    % assume already cellstr-ish
    c = v(:);
    % coerce entries to char
    for i=1:numel(c)
        c{i} = char(string(c{i}));
    end
else
    c = cellstr(string(v(:)));
end
if nargin > 1 && numel(c) ~= n
    % try broadcasting single value
    if numel(c)==1
        c = repmat(c, n, 1);
    else
        % pad/truncate
        c = c(1:min(end,n));
        if numel(c) < n
            c(end+1:n,1) = {''};
        end
    end
end
end


function s = mode_string(x)
% mode for strings/cellstr/numeric -> string
if isstring(x)
    u = unique(x);
    if isempty(u), s=""; else, s=u(1); end
elseif iscell(x)
    u = unique(string(x));
    if isempty(u), s=""; else, s=u(1); end
else
    x = x(:);
    x = x(isfinite(x));
    if isempty(x), s=""; else, s=string(mode(x)); end
end
end



% ============================================================
% ===================== MATH UTILITIES =======================
% ============================================================

function s = stripChanSuffix(fname)
% Remove final _NN.mat
% Example: abc_def_03.mat -> abc_def
s = regexprep(string(fname), "_\d\d\.mat$", "");
end


function out = sanitize_for_filename(s)
s = string(s);
out = regexprep(s, '[^\w\-]+', '_');     % keep letters/numbers/_/-
out = regexprep(out, '_+', '_');
out = strip(out, "_");
if strlength(out)==0, out="NA"; end
end


function m = robust_mad(x)
x = x(:);
x = x(isfinite(x));
if isempty(x), m = NaN; return; end
med = median(x);
m = median(abs(x-med));
end


function z = robust_z(x)
x = x(:);
med = median(x,'omitnan');
m = robust_mad(x);
if ~isfinite(m) || m==0
    z = (x - med);
else
    z = (x - med) ./ (1.4826*m);
end
end


function [xx, dd] = simple_kde(x, nGrid)
% Simple Gaussian KDE without toolboxes
if length(nGrid) > 1
    nGrid = length(nGrid);
end
x = x(:);
x = x(isfinite(x));
if numel(x) < 3
    xx = linspace(min(x), max(x), nGrid);
    dd = zeros(size(xx));
    return
end

xmin = prctile(x, 1);
xmax = prctile(x, 99);
if xmin==xmax
    xmin = min(x); xmax = max(x) + eps;
end
xx = linspace(xmin, xmax, nGrid);

% Silverman's rule bandwidth
sx = std(x,0);
if ~isfinite(sx) || sx==0, sx = robust_mad(x)/0.6745; end
h = 1.06 * sx * numel(x)^(-1/5);
if ~isfinite(h) || h<=0
    h = (xmax-xmin)/30;
end

% KDE
dd = zeros(size(xx));
c = 1/(sqrt(2*pi)*h*numel(x));
for i=1:numel(x)
    dd = dd + exp(-0.5*((xx - x(i))./h).^2);
end
dd = dd * c;
end


function s = chan_label(chanDat)
try
    s = string(chanDat.labels{chanDat.chi});
catch
    s = "chan"+string(chanDat.chi);
end
end


function y = ternary(cond, a, b)
if cond, y=a; else, y=b; end
end
