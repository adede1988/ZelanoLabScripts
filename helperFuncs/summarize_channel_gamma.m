
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

