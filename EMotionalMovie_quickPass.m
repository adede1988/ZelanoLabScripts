function [] = EMotionalMovie_quickPass()
codePre = 'G:\My Drive\GitHub\';
datPre  = {'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ...
         'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\', ...
         'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

sessionIDs = {'251009_OBE_NWU_CP_1', ... 
           '250225_OBE_NWU_AS_4', ...  
            '250904_OBE_NWU_TI'};           
datPrei = [2,2,2];       

addpath(genpath([codePre 'ZelanoLabScripts']))
addpath(genpath([codePre 'myFrequentUse']))

for s = 2:numel(sessionIDs)
    S.id   = sessionIDs{s};
    S.root = datPre{datPrei(s)};

    outDat = load(fullfile(S.root, S.id, 'preProc', ...
                [S.id '_' 'EmotionalMovieTask' 'preproc.mat']));
    if isfield(outDat, 'outDat')
        outDat = outDat.outDat;
    end

    idx = cellfun(@(x) contains(x, 'macBP'), outDat.labels);
    macOut = outDat.data(idx,:); 
    %look at channel 3 for all for now: 
    chanOfInterest = 3; 
    %gamma frequencies for search: 
    gamFrex = logspace(log10(25), log10(70), 100);
    %pull out respiration data: 
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx,:); 
    rspDat = rspDat(outDat.rspIDX,:);
    rspDat = rspDat .* outDat.rspFlip;
    %regenerate the typeidx vector for condition binning: 
    typeidx = zeros(size(outDat.data,2),1);
    tim = [1/outDat.fs:1/outDat.fs:size(outDat.data,2)*(1/outDat.fs)];
    for ii = 1:length(outDat.TTL.timStamp)-1
        outDat.TTL.timStampEnd(ii) = outDat.TTL.timStamp(ii+1) - ...
                                            outDat.fs*(1.6);
        typeidx(outDat.TTL.timStamp(ii):outDat.TTL.timStamp(ii+1)) =...
                        outDat.TTL.type(ii); 
    end
    typeidx = typeidx'; 

    %% =========================== USER INPUTS ===============================
    figureSavePath = fullfile(outDat.figs);  % change if desired
    if ~exist(figureSavePath, 'dir'); mkdir(figureSavePath); end
    
    chanOfInterest = 3;          % your best anatomical channel (fixed)
    gamBand = [25 70];           % Hz
    gamFrex = linspace(gamBand(1), gamBand(2), 181);  % fine grid for peak search
    burstPct = 95;               % percentile threshold
    minCycles = 2;               % min cycles for gamma burst
    maxBurst_s = 1.0;            % max allowed burst length (s)
    trialBuff_s = 2;             % trial window on each side of onset
    tfFrex = logspace(log10(2), log10(150), 100);   % TF freqs (2–150 Hz)
    phaseAtOnsetFreqs = [2 5 8 12 15 20];           % Hz (F, G)
    phaseAfterKcycles = [0 1 2 3 4 5 8 12 15 20];           % cycles after onset (C)
    hiFreqAgg = [50 150];        % Hz band for aggregate z-power (H)
    semAlpha = 0.25;             % SEM shading opacity
    nPhaseBins = 32;             % for polar “frequency polygons”
    nLenBins = 40;               % for burst length distributions
    % rng(1);                      % for any stochastic cleaning reproducibility (if used)
    
    % Inputs expected: rspDat (1 x time), macOut (channels x time), typeidx (1 x time)
    %                   outDat.fs (Hz), outDat.labels (cellstr)
    fs = outDat.fs;
    
    %% ===================== 0) BASIC SIGNALS & GUARDS =======================
    x = double(squeeze(macOut(chanOfInterest, :)));
    T = numel(x);
    assert(isvector(rspDat) && numel(rspDat)==T, 'rspDat must be 1 x time matching ephys.');
    assert(isvector(typeidx) && numel(typeidx)==T, 'typeidx must be 1 x time matching ephys.');
    typeidx = double(typeidx(:)).';  % ensure row vector
    
    % Optional: soft detrend / HPF could go here if needed (left as-is to preserve your pipeline)
    
    %% ================= 1) GAMMA PEAK, BURSTS & RESP PHASE ==================
    % Time-frequency at gamma grid (whole recording, single channel)
    [phase_gam, pow_gam] = multiphasevec3(gamFrex, x, fs, 4, 0);
    % pow_gam assumed dims: (1, nFrex, nTime) or similar; use squeeze robustly
    pow_gam = squeeze(pow_gam);   % => nFrex x nTime or nTime x nFrex; normalize shape
    if size(pow_gam,1) ~= numel(gamFrex)
        pow_gam = permute(pow_gam, [2 1]); % make it (nFrex x nTime)
    end
    phase_gam = squeeze(phase_gam); % make it (1|n?, nFrex, nTime)
    if ndims(phase_gam)==3
        % retain as (nFrex x nTime)
        phase_gam = squeeze(phase_gam(1,:,:));
    end
    if size(phase_gam,1) ~= numel(gamFrex)
        phase_gam = permute(phase_gam, [2 1]); % (nFrex x nTime)
    end
    
    % Power spectrum (A): average across time
    meanPowSpec = mean(pow_gam, 2);                          % (nFrex x 1)
    [~, maxidx]  = max(meanPowSpec);
    peakGamHz    = gamFrex(maxidx);
    Tgam         = 1/peakGamHz;
    
    % Figure A: gamma power spectrum with annotation
    fhA = figure('Color','w', 'position', [0,0, 500, 500], 'visible', false); hold on
    plot(gamFrex, meanPowSpec, 'LineWidth', 2);
    xline(peakGamHz, 'r--', 'LineWidth', 2);
    text(peakGamHz, max(meanPowSpec)*0.95, sprintf('  Peak = %.2f Hz', peakGamHz), ...
        'Color','r','FontWeight','bold','VerticalAlignment','top');
    xlabel('Frequency (Hz)'); ylabel('Power (a.u.)');
    title(sprintf('%s: Channel %d Gamma-band Spectrum', ...
        outDat.sessID, chanOfInterest), 'interpreter', 'none');
    tightfig();
    saveas(fhA, fullfile(figureSavePath, sprintf('A_gamma_spectrum_chan%d.jpg', chanOfInterest)));
    
    % Burst detection at peak frequency
    gamPowPeak = pow_gam(maxidx, :);                 % (1 x time)
    gamPowPeak_sm = smoothdata(gamPowPeak, 'gaussian', round(0.1*fs));  % ~100 ms
    burstThresh = prctile(gamPowPeak_sm, burstPct);
    
    % Upward crossings (onsets)
    onsets = find(gamPowPeak_sm(1:end-1) < burstThresh & gamPowPeak_sm(2:end) >= burstThresh) + 1;
    
    % Align onset to nearest gamma phase peak within ±25 samples
    pad = 25;
    onsets = onsets(onsets > pad & onsets <= T-pad);
    phPeak = phase_gam(maxidx, :);  % phase at peak gamma freq over time
    alignFun = @(ix) ix - pad + argmin_phase_distance(phPeak(ix-pad:ix+pad), pi/2); % nearest +pi/2
    onsets = arrayfun(alignFun, onsets);
    
    % Remove first/last few to avoid edges on later trialing
    onsets = onsets(onsets > trialBuff_s*fs & onsets <= (T - trialBuff_s*fs));
    
    % Drop bursts in bad epochs (if present)
    badLabIdx = find(cellfun(@(s) strcmpi(s,'badTS'), outDat.labels));
    if ~isempty(badLabIdx)
        badTS = logical(outDat.data(badLabIdx, onsets));
        onsets = onsets(~badTS);
    end
    
    % Compute burst lengths (until power drops below threshold again)
    offIdx = arrayfun(@(ix) ix - 1 + find(gamPowPeak_sm(ix:end-1) >= burstThresh & ...
                                          gamPowPeak_sm(ix+1:end) <  burstThresh, 1, 'first'), onsets);
    % sanitize (if trailing no-off found)
    valid = ~cellfun(@isempty, arrayfun(@(a,b) {b}, onsets, offIdx));
    onsets = onsets(valid);
    offIdx = offIdx(valid);
    burstLen_s = (offIdx - onsets) / fs;
    
    % Length constraints: >= 2 cycles of peak gamma and <= 1 s
    minLen_s = minCycles * Tgam;
    keep = burstLen_s >= minLen_s & burstLen_s <= maxBurst_s;
    onsets = onsets(keep);
    burstLen_s = burstLen_s(keep);
    
    % Behavioral types per burst
    burstType = typeidx(onsets);
    
    % Respiration phase at burst onset (E)
    rspPhase = angle(hilbert(smoothdata(double(rspDat(:)).', 'gaussian', fs)));
    burstRspPhase = rspPhase(onsets);
    
    %% ============ 2) TRIALIZATION, CLEANING, ERP PER CONDITION ============
    W = round(trialBuff_s * fs);
    nBursts = numel(onsets);
    nChan = size(macOut,1);
    trialDat = zeros(nChan, 2*W, nBursts, 'double');
    trialRsp = zeros(2*W, nBursts, 'double');
    
    for k = 1:nBursts
        t0 = onsets(k);
        trialDat(:,:,k) = macOut(:, t0-W+1 : t0+W);  % inclusive 2W samples
        trialRsp(:,k)   = rspDat(  t0-W+1 : t0+W);
    end
    
    % Covariance-matrix cleaning (you already have this)
    badIDX = covMatClean_breathing(trialDat, ones(nBursts,1), ones(nBursts,1) * (2*W));
    keepTrial = true(1,nBursts);
    keepTrial(badIDX) = false;
    
    trialDat = trialDat(:,:,keepTrial);
    trialRsp = trialRsp(:,keepTrial);
    onsets   = onsets(keepTrial);
    burstLen_s = burstLen_s(keepTrial);
    burstType  = burstType(keepTrial);
    burstRspPhase = burstRspPhase(keepTrial);
    nBursts = sum(keepTrial);
    
    % ERP per condition (B) on channel of interest
    tRel = ((-W+1):(W)) / fs;  % seconds relative to onset
    Cnames = {'Neutral','Happy','Sad'};
    Ccols  = lines(3);  % distinct colors
    erp = cell(1,3); sem = cell(1,3); nC = zeros(1,3);
    for c = 1:3
        idxC = find(burstType==c);
        nC(c) = numel(idxC);
        Y = squeeze(trialDat(chanOfInterest, :, idxC)).'; % trials x time
        erp{c} = mean(Y, 1);
        sem{c} = std(Y, [], 1) ./ sqrt(max(1, size(Y,1)));
    end
    
    fhB = figure('Color','w', 'position', [0,0,1200,600], 'visible', false); hold on
    for c = 1:3
        plot(tRel, erp{c}, 'Color', Ccols(c,:), 'LineWidth', 2);
    end
    xline(0,'k--'); grid on
    xlabel('Time from gamma-onset (s)'); ylabel('Voltage (µV)');
    title(sprintf('%s: ERP (chan %d), cleaned trials',outDat.sessID, ...
                chanOfInterest), 'interpreter', 'none');
    legend(sprintf('%s (n=%d)',Cnames{1},nC(1)), ...
           sprintf('%s (n=%d)',Cnames{2},nC(2)), ...
           sprintf('%s (n=%d)',Cnames{3},nC(3)), 'Location','best',...
           'AutoUpdate','off');
    for c = 1:3
        shaded_sem(tRel, erp{c}, sem{c}, Ccols(c,:), semAlpha);
    end
    tightfig();
    saveas(fhB, fullfile(figureSavePath, sprintf('B_ERP_sem_chan%d.jpg', chanOfInterest)));
    
    %% ======= 3) PHASE OF GAMMA AFTER K CYCLES FROM ONSET (C; polar) =======
    % Use gamma phase series at peak frequency (computed earlier as phPeak)
    % For each burst: sample phase at onset + K*Tgam (nearest sample)
    make_phase_polygons_after_cycles(phPeak, onsets, burstType, Tgam, phaseAfterKcycles, ...
        Cnames, figureSavePath, 'C_gammaPhase_afterKcycles', outDat.fs);
    
    %% ======== 4) BURST LENGTH DISTRIBUTION PER CONDITION (D; linear) =======
    fhD = figure('Color','w'); hold on
    for c = 1:3
        vals = burstLen_s(burstType==c);
        [ctrs, dens] = freqpoly(vals, nLenBins);
        plot(ctrs, dens, 'LineWidth', 2, 'Color', Ccols(c,:));
    end
    xlabel('Gamma burst length (s)'); ylabel('Density');
    title('Gamma burst length distribution by condition');
    legend(Cnames, 'Location','best'); grid on
    tightfig();
    saveas(fhD, fullfile(figureSavePath, 'D_burst_length_distribution.jpg'));
    
    %% === 5) RESPIRATION PHASE AT BURST ONSET PER CONDITION (E; polar) ======
    make_phase_polygons_polar(burstRspPhase, burstType, Cnames, ...
        figureSavePath, 'E_resp_phase_at_gamma_onset', nPhaseBins, 'Respiration Phase at Gamma Burst');
    
    %% === 6) OTHER-FREQ PHASES AT GAMMA ONSET: 2,5,8,12,15,20 Hz (F; polar)
    [phase_sel, ~] = multiphasevec3(phaseAtOnsetFreqs, x, fs, 4, 0);
    phase_sel = squeeze(phase_sel); % expect (nFreq x nTime)
    if size(phase_sel,1) ~= numel(phaseAtOnsetFreqs)
        phase_sel = permute(phase_sel, [2 1]);
    end
    for fi = 1:numel(phaseAtOnsetFreqs)
        thisPhase = phase_sel(fi, onsets);
        myTitle = ['phase distribution for ' num2str(phaseAtOnsetFreqs(fi)) ' Hz at time of gamma peak']; 
        make_phase_polygons_polar(thisPhase, burstType, Cnames, ...
            figureSavePath, sprintf('F_phase_%gHz_at_onset', phaseAtOnsetFreqs(fi)), nPhaseBins, myTitle);
    end
    
    %% === 7) Z-SCORED POWER DISTRIBUTIONS @ (2,5,8,12,15,20 Hz) (G; linear)
    [~, pow_sel] = multiphasevec3(phaseAtOnsetFreqs, x, fs, 4, 0);
    pow_sel = squeeze(pow_sel); % (nFreq x nTime)
    if size(pow_sel,1) ~= numel(phaseAtOnsetFreqs)
        pow_sel = permute(pow_sel, [2 1]);
    end
    % z-score within frequency across full time
    parfor fi = 1:numel(phaseAtOnsetFreqs)
        pow_sel(fi,:) = myChanZscore(pow_sel(fi,:));
    end
    for fi = 1:numel(phaseAtOnsetFreqs)
        vals = pow_sel(fi, onsets);  % z-power at onset
        fhG = figure('Color','w'); hold on
        for c = 1:3
            v = vals(burstType==c);
            [ctrs, dens] = freqpoly(v, round(sqrt(numel(v))+10));
            plot(ctrs, dens, 'LineWidth', 2, 'Color', Ccols(c,:));
        end
        xlabel('Z-scored power'); ylabel('Density'); grid on
        title(sprintf('Z-power distribution @ %g Hz (by condition)', phaseAtOnsetFreqs(fi)));
        legend(Cnames, 'Location','best');
        tightfig();
        saveas(fhG, fullfile(figureSavePath, sprintf('G_zpower_%gHz.jpg', phaseAtOnsetFreqs(fi))));
    end
    
    %% === 8) HI-FREQ (50–150 Hz) MEAN Z-POWER in ±0.5 s AROUND ONSET (H) ===
    hiIdx = tfFrex >= hiFreqAgg(1) & tfFrex <= hiFreqAgg(2);
    % Compute TF on whole signal for speed; then window per burst
    tmp = squeeze(trialDat(chanOfInterest,:,:)); 
    [phase_whole, pow_whole] = multiphasevec3(tfFrex, tmp', fs, 6, 0); % better cycles for hi-freq
    % pow_whole = squeeze(pow_whole);  % (nFreq x nTime)
    % if size(pow_whole,1) ~= numel(tfFrex)
    %     pow_whole = permute(pow_whole, [2 1]);
    % end
    
    % z-score within freq
    parfor fi = 1:numel(tfFrex)
     fi
     tic
        tmp = squeeze(pow_whole(:,fi,:))';
        tmp = tmp(fs:fs*3,:); 
        tmp = myChanZscore(tmp);
        slice = squeeze(pow_whole(:,fi,:))';
        slice(fs:fs*3,:) = tmp; 
        pow_whole(:,fi,:) = slice';
        toc
    end
    halfWin = round(0.5 * fs);
    valsH = mean(pow_whole(:,hiIdx, ...
                round((trialBuff_s-.5)*fs) : round((trialBuff_s+.5)*fs) ),...
                                [2,3]);
   
    
    fhH = figure('Color','w'); hold on
    for c = 1:3
        v = valsH(burstType==c);
        [ctrs, dens] = freqpoly(v, round(sqrt(numel(v))+10));
        plot(ctrs, dens, 'LineWidth', 2, 'Color', Ccols(c,:));
    end
    xlabel(sprintf('Mean z-power %d–%d Hz in ±0.5 s', hiFreqAgg(1), hiFreqAgg(2)));
    ylabel('Density'); grid on; legend(Cnames, 'Location','best');
    title('High-frequency aggregate power around gamma onset (by condition)');
    tightfig();
    saveas(fhH, fullfile(figureSavePath, 'H_hiFreqAgg_pm500ms.jpg'));
    
    

    %% === 9) Count of gamma events per condition type (I) ===
    conds = 1:3;
    counts = arrayfun(@(c) sum(burstType==c), conds);
    Ntot   = sum(counts);
    perc   = 100 * counts / max(1,Ntot);

    % Plot
    fhI = figure('Color','w'); 
    b = bar(conds, counts, 'FaceColor','flat', 'EdgeColor','none', 'BarWidth',0.6);
    for i = 1:3
        b.CData(i,:) = Ccols(i,:);
    end
    xticks(conds); xticklabels(Cnames);
    ylabel('Number of gamma events');
    title('Gamma events per condition');
    grid on

    % Add count + percentage labels above bars
    yl = ylim;
    for i = 1:3
        text(i, counts(i) + 0.02*range(yl), sprintf('%d (%.1f%%)', counts(i), perc(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
    end
    ylim([0, max(counts)*1.15 + eps]); % a little headroom for labels

    % Save
    tightfig(); 
    saveas(fhI, fullfile(figureSavePath, 'I_gamma_event_counts.jpg'));


%% === 10) Power time frequency plots for each condition (J) === 
halfWin_s   = 1;                      % ±1 second
halfWin_smp = round(halfWin_s * fs);  % samples
centerIdx   = round((size(pow_whole,3) + 1) / 2); % assume trials centered on onset

startIdx = centerIdx - halfWin_smp;
endIdx   = centerIdx + halfWin_smp;



tIdx = startIdx:endIdx;
tRel = (tIdx - centerIdx) / fs;       % time (s) relative to onset
nTimeShort = numel(tIdx);

% Slice into shorter window: nEvents x nFreq x nTimeShort
pow_short = pow_whole(:, :, tIdx);

% ------------ Condition means & global color limits ------------
condMeans = cell(1,3);
clims = [inf -inf];

for c = 1:3
    idxC = find(burstType == c);
    if isempty(idxC)
        condMeans{c} = nan(nFreq, nTimeShort);
        warning('No events for condition %d (%s).', c, Cnames{c});
    else
        % Average over events: result nFreq x nTimeShort
        condMeans{c} = squeeze(mean(pow_short(idxC, :, :), 1, 'omitnan'));
    end
    cm = condMeans{c};
    clims(1) = -5;
    clims(2) = 5;
end

% Fallback color limits if all-NaN or flat
if ~all(isfinite(clims)) || diff(clims)==0
    clims = [min(-1, clims(1)) max(1, clims(2))];
end

% ------------ Individual condition figures & saving ------------
for c = 1:3
    fh = figure('Color','w','Position',[150 150 540 460]);
    imagesc(tRel, [], condMeans{c});
    yticks([1:10:100])
    yticklabels(round(tfFrex(1:10:100)))
    set(gca, 'YDir', 'normal');
    caxis(clims);
    colormap(parula);
    title(sprintf('Time–frequency power: %s (±1 s)', Cnames{c}));
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    xline(0,'k--','LineWidth',1.25);
    cb = colorbar;
    cb.Label.String = 'Power (a.u.)';
    try, tightfig(); catch, end

    fname = sprintf('TF_power_%s_pm1s.jpg', lower(Cnames{c}));
    saveas(fh, fullfile(figureSavePath, fname));
end

%% === 11) ITPC time frequency plots for each condition (K) === 
halfWin_s   = 1;                      % ±1 second
halfWin_smp = round(halfWin_s * fs);  % samples
centerIdx   = round((size(pow_whole,3) + 1) / 2); % assume trials centered on onset

startIdx = centerIdx - halfWin_smp;
endIdx   = centerIdx + halfWin_smp;



tIdx = startIdx:endIdx;
tRel = (tIdx - centerIdx) / fs;       % time (s) relative to onset
nTimeShort = numel(tIdx);

% Slice into shorter window: nEvents x nFreq x nTimeShort
phase_short = phase_whole(:, :, tIdx);

% ------------ Condition means & global color limits ------------
condMeans = cell(1,3);
clims = [inf -inf];

for c = 1:3
    idxC = find(burstType == c);
    if isempty(idxC)
        condMeans{c} = nan(nFreq, nTimeShort);
        warning('No events for condition %d (%s).', c, Cnames{c});
    else
        % Average over events: result nFreq x nTimeShort
        condMeans{c} = squeeze(abs(mean(exp(1i*phase_short(idxC, :, :)), 1)));
    end
    cm = condMeans{c};
    clims(1) = 0;
    clims(2) = .1;
end

% Fallback color limits if all-NaN or flat
if ~all(isfinite(clims)) || diff(clims)==0
    clims = [min(-1, clims(1)) max(1, clims(2))];
end

% ------------ Individual condition figures & saving ------------
for c = 1:3
    fh = figure('Color','w','Position',[150 150 540 460]);
    imagesc(tRel, [], condMeans{c});
    yticks([1:10:100])
    yticklabels(round(tfFrex(1:10:100)))
    set(gca, 'YDir', 'normal');
    caxis(clims);
    colormap(parula);
    title(sprintf('Time–frequency ITPC: %s (±1 s)', Cnames{c}));
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    xline(0,'k--','LineWidth',1.25);
    cb = colorbar;
    cb.Label.String = 'Power (a.u.)';
    try, tightfig(); catch, end

    fname = sprintf('TF_ITPC_%s_pm1s.jpg', lower(Cnames{c}));
    saveas(fh, fullfile(figureSavePath, fname));
end



%% ==============================================================
% EXTRA FIGURES FOR GAMMA–RESPIRATION PIPELINE
% ==============================================================

fs   = outDat.fs;
x    = double(squeeze(macOut(chanOfInterest, :)));
T    = numel(x);
Cnames = {'Neutral','Happy','Sad'};
Ccols  = lines(3);

if ~exist(figureSavePath,'dir'); mkdir(figureSavePath); end

%% 1a) Six randomly chosen raw & gamma-band filtered traces (2 per condition)
% -------------------------------------------------------------------------
                       % reproducible
winRaw_s   = 1;                % ±1 s
winRaw_smp = round(winRaw_s*fs);

for c = 1:3
    idxC = find(burstType==c);
    if isempty(idxC), continue; end
    nToPlot = min(10, numel(idxC));
    idxSel = idxC(randperm(numel(idxC), nToPlot));

    for k = 1:nToPlot
        bIdx = idxSel(k);
        t0   = onsets(bIdx);
        if t0 <= winRaw_smp || t0 > (T - winRaw_smp), continue; end

        segIdx = (t0-winRaw_smp) : (t0+winRaw_smp);
        tRel   = (segIdx - t0) / fs;

        rawSeg   = x(segIdx);
        gamSeg   = bandpass(rawSeg,[25 70],fs);  % gamma-band filtered

        fh = figure('Color','w','Position',[100 100 600 400]);
        
        plot(tRel, rawSeg, 'k'); hold on
        xline(0,'r--','LineWidth',1.5);
        xlabel('Time (s)'); ylabel('Voltage (µV)');
        % title(sprintf('%s: raw (event %d)', Cnames{c}, bIdx));

        
        plot(tRel, gamSeg, 'Color',[0.2 0.4 0.8]); hold on
        xline(0,'r--','LineWidth',1.5);
        xlabel('Time (s)'); ylabel('Gamma-filtered (a.u.)');
        % title(sprintf('%s: gamma-band (25–70 Hz)', Cnames{c}));

        sgtitle(sprintf('Burst-centered trace (chan %d)', chanOfInterest));

        fname = sprintf('rawSingleEventTrace_%s_evt%04d.jpg', lower(Cnames{c}), bIdx);
        saveas(fh, fullfile(figureSavePath, fname));
    end
end

%% 2a) Freqpoly of inter-burst intervals (IBI) per condition
% -------------------------------------------------------------------------
fh = figure('Color','w'); hold on
edges = 0:0.1:5;              % 0–5 s, 100 ms bins

for c = 1:3
    ts  = sort(onsets(burstType==c) / fs);  % in seconds
    if numel(ts) < 2, continue; end
    ibi = diff(ts);                          % inter-burst intervals

    [cts, e] = histcounts(ibi, edges, 'Normalization','pdf');
    ctrs     = e(1:end-1) + diff(e)/2;
    plot(ctrs, cts, 'LineWidth',2, 'Color', Ccols(c,:));
end
xlabel('Inter-burst interval (s)');
ylabel('Density');
title('Inter-burst intervals by condition');
legend(Cnames, 'Location','best');
grid on
saveas(fh, fullfile(figureSavePath,'extra_IBI_freqpoly_by_condition.jpg'));


%% 2a – autocorrelogram of burst times (±5 s) per condition
% -------------------------------------------------------------------------
lagMax_s = 5;
lagEdges = -lagMax_s:0.1:lagMax_s;

fh = figure('Color','w','Position',[100 100 1000 350]);
hold on 
for c = 1:3
    ts = sort(onsets(burstType==c) / fs);   % in seconds
    allDiffs = [];

    for i = 1:numel(ts)
        d = ts - ts(i);
        d = d(d ~= 0 & abs(d) <= lagMax_s);
        allDiffs = [allDiffs; d(:)];
    end

    if ~isempty(allDiffs)
        [cts, e] = histcounts(allDiffs, lagEdges, 'Normalization','pdf');
        ctrs = e(1:end-1) + diff(e)/2;
        plot(ctrs, cts, 'LineWidth',2, 'Color', Ccols(c,:));
        hold on; 
    end
    xlabel('Lag (s)');
    ylabel('Density');
    
    grid on
end
legend(Cnames, 'AutoUpdate','off')
xline(0,'k--','LineWidth',1);
saveas(fh, fullfile(figureSavePath,'extra_burst_autocorr_pm5s.jpg'));

%% 2a – mean theta (4–8 Hz) power around gamma peaks (±2 s) per condition
% -------------------------------------------------------------------------

bandNames = {'delta', 'theta', 'alpha', 'beta'}; 
bandVals = {[2 4], [4 8], [8 15], [15 25]};

for bi = 1: 4
    thetaSig = bandpass(x,bandVals{bi},fs);
    thetaPow = abs(hilbert(thetaSig)).^2;
    
    winTheta_s   = 2;
    winTheta_smp = round(winTheta_s * fs);
    tTheta       = (-winTheta_smp:winTheta_smp) / fs;
    Ltheta       = numel(tTheta);
    
    thetaMat = nan(numel(onsets), Ltheta);
    
    validBurst = onsets > winTheta_smp & onsets <= (T - winTheta_smp);
    for k = find(validBurst(:))'
        t0 = onsets(k);
        idx = (t0-winTheta_smp):(t0+winTheta_smp);
        thetaMat(k,:) = thetaPow(idx);
    end
    
    fh = figure('Color','w'); hold on
    for c = 1:3
        idxC = find(burstType==c & validBurst(:)');
        if isempty(idxC), continue; end
        Y = thetaMat(idxC,:);
        m = mean(Y,1,'omitnan');
        se = std(Y, [], 1,'omitnan') ./ sqrt(size(Y,1));
        % shaded mean
        
        plot(tTheta, m, 'Color',Ccols(c,:), 'LineWidth',2);
    end
    xline(0,'k--','LineWidth',1.5);
    xlabel('Time from gamma onset (s)');
    ylabel([bandNames{bi} 'power (a.u.)']);
    title([bandNames{bi} ' power around gamma bursts']);
    legend(Cnames, 'Location','best', 'autoupdate', false);
    for c = 1:3
        idxC = find(burstType==c & validBurst(:)');
        if isempty(idxC), continue; end
        Y = thetaMat(idxC,:);
        m = mean(Y,1,'omitnan');
        se = std(Y, [], 1,'omitnan') ./ sqrt(size(Y,1));
        patch([tTheta fliplr(tTheta)], [m+se fliplr(m-se)], Ccols(c,:), ...
              'FaceAlpha',0.2,'EdgeColor','none');
    end
    grid on
    saveas(fh, fullfile(figureSavePath,['extra_' bandNames{bi} '_power_around_gamma_pm2s.jpg']));
end


%% 3a) Burst-onset respiration phase vs burst length (scatter)
% -------------------------------------------------------------------------
% burstRspPhase: phase at onset (radians); burstLen_s: burst length (s)

fh = figure('Color','w'); hold on
for c = 1:3
    idxC = burstType==c;
    scatter(rad2deg(burstRspPhase(idxC)), burstLen_s(idxC), ...
        15, Ccols(c,:), 'filled', 'MarkerFaceAlpha',0.4);
end
xlabel('Respiration phase at gamma onset (degrees)');
ylabel('Gamma burst length (s)');
title('Burst length vs respiration phase at onset');
legend(Cnames, 'Location','best');
xlim([-180 180]);
grid on
saveas(fh, fullfile(figureSavePath,'extra_burstLen_vs_respPhase_scatter.jpg'));


%% 3a – mean gamma power locked to respiration onset (-2 to +5 s)
% -------------------------------------------------------------------------
respOnsets = outDat.behDat.finalOnset(:);      % indices into full time
respCond   = outDat.behDat.condition(:);       % 1/2/3 per respiration

% Recompute gamma power time series at peakGamHz if needed
[~, powGam1] = multiphasevec3(peakGamHz, x, fs, 6, 0);
powGam1 = squeeze(powGam1);          % 1 x T or T x 1
if numel(powGam1) ~= T
    powGam1 = powGam1(:).';
end
gamPow_ts = powGam1;                 % gamma power trace

winPre_s   = 2;
winPost_s  = 5;
preSmp     = round(winPre_s*fs);
postSmp    = round(winPost_s*fs);
tResp      = (-preSmp:postSmp) / fs;
Lresp      = numel(tResp);

gamRespMat = nan(numel(respOnsets), Lresp);

validResp = respOnsets > preSmp & respOnsets <= (T - postSmp);
for i = find(validResp(:))'
    t0  = respOnsets(i);
    idx = (t0-preSmp):(t0+postSmp);
    gamRespMat(i,:) = gamPow_ts(idx);
end

fh = figure('Color','w'); hold on
for c = 1:3
    idxC = find(respCond==c & validResp);
    if isempty(idxC), continue; end
    Y = gamRespMat(idxC,:);
    m = mean(Y,1,'omitnan');
    se = std(Y,[],1,'omitnan') ./ sqrt(size(Y,1));
   
    plot(tResp, m, 'Color',Ccols(c,:), 'LineWidth',2);
end
xline(0,'k--','LineWidth',1.5);
xlabel('Time from respiration onset (s)');
ylabel(sprintf('Gamma power @ %.1f Hz (a.u.)', peakGamHz));
title('Gamma power locked to respiration onset');
legend(Cnames, 'Location','best', 'autoupdate', false);
grid on
for c = 1:3
    idxC = find(respCond==c & validResp);
    if isempty(idxC), continue; end
    Y = gamRespMat(idxC,:);
    m = mean(Y,1,'omitnan');
    se = std(Y,[],1,'omitnan') ./ sqrt(size(Y,1));
    patch([tResp fliplr(tResp)], [m+se fliplr(m-se)], Ccols(c,:), ...
          'FaceAlpha',0.2,'EdgeColor','none');
end
saveas(fh, fullfile(figureSavePath,'extra_gamma_power_locked_to_resp_onset.jpg'));

%% 4a) Time–frequency difference maps (Happy–Neutral, Sad–Neutral)
%     using pow_whole (events x freq x time) and ±1 s window
% -------------------------------------------------------------------------
[nEvents, nFreq, nTime] = size(pow_whole);
centerIdx   = round((nTime + 1)/2);
halfWin_s   = 1;
halfWin_smp = round(halfWin_s * fs);
startIdx    = max(1, centerIdx-halfWin_smp);
endIdx      = min(nTime, centerIdx+halfWin_smp);
tIdx        = startIdx:endIdx;
tRelTF      = (tIdx - centerIdx) / fs;
Ltf         = numel(tIdx);

pow_short = pow_whole(:, :, tIdx);   % events x freq x timeShort

condMeans = cell(1,3);
for c = 1:3
    idxC = find(burstType==c);
    if isempty(idxC)
        condMeans{c} = nan(nFreq, Ltf);
    else
        condMeans{c} = squeeze(mean(pow_short(idxC,:,:),1,'omitnan'));
    end
end

diff_HN = condMeans{2} - condMeans{1};   % Happy - Neutral
diff_SN = condMeans{3} - condMeans{1};   % Sad   - Neutral

% Symmetric color limits across both maps
mx = max(abs([diff_HN(:); diff_SN(:)]), [], 'omitnan');
if ~isfinite(mx) || mx==0, mx = 1; end
cl = [-mx mx];

% Happy - Neutral
fh = figure('Color','w','Position',[150 150 540 460]);
imagesc(tRelTF, [], diff_HN);
yticks([1:10:100])
    yticklabels(round(tfFrex(1:10:100)))
set(gca,'YDir','normal'); caxis([-5, 10]); colormap(parula);  % requires redblue colormap or change to parula
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('TF difference: Happy – Neutral');
xline(0,'k--','LineWidth',1.25);
cb = colorbar; cb.Label.String = 'Power diff (a.u.)';
saveas(fh, fullfile(figureSavePath,'extra_TF_diff_Happy_minus_Neutral_pm1s.jpg'));

% Sad - Neutral
fh = figure('Color','w','Position',[150 150 540 460]);
imagesc(tRelTF, [], diff_SN);
yticks([1:10:100])
    yticklabels(round(tfFrex(1:10:100)))
set(gca,'YDir','normal'); caxis([-5, 10]); 
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('TF difference: Sad – Neutral');
xline(0,'k--','LineWidth',1.25);
cb = colorbar; cb.Label.String = 'Power diff (a.u.)';
saveas(fh, fullfile(figureSavePath,'extra_TF_diff_Sad_minus_Neutral_pm1s.jpg'));


%% 4a) Time–frequency ITPC difference maps (Happy–Neutral, Sad–Neutral)
%     using pow_whole (events x freq x time) and ±1 s window
% -------------------------------------------------------------------------
[nEvents, nFreq, nTime] = size(phase_whole);
centerIdx   = round((nTime + 1)/2);
halfWin_s   = 1;
halfWin_smp = round(halfWin_s * fs);
startIdx    = max(1, centerIdx-halfWin_smp);
endIdx      = min(nTime, centerIdx+halfWin_smp);
tIdx        = startIdx:endIdx;
tRelTF      = (tIdx - centerIdx) / fs;
Ltf         = numel(tIdx);

pow_short = phase_whole(:, :, tIdx);   % events x freq x timeShort

condMeans = cell(1,3);
for c = 1:3
    idxC = find(burstType==c);
    if isempty(idxC)
        condMeans{c} = nan(nFreq, Ltf);
    else
        condMeans{c} = squeeze(abs(mean(exp(1i*pow_short(idxC,:,:)),1) ) );
    end
end

diff_HN = condMeans{2} - condMeans{1};   % Happy - Neutral
diff_SN = condMeans{3} - condMeans{1};   % Sad   - Neutral

% Symmetric color limits across both maps
mx = max(abs([diff_HN(:); diff_SN(:)]), [], 'omitnan');
if ~isfinite(mx) || mx==0, mx = 1; end
cl = [-mx mx];

% Happy - Neutral
fh = figure('Color','w','Position',[150 150 540 460]);
imagesc(tRelTF, [], diff_HN);
yticks([1:10:100])
    yticklabels(round(tfFrex(1:10:100)))
set(gca,'YDir','normal'); caxis([-.1 .1]); colormap(parula);  % requires redblue colormap or change to parula
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('TF difference: Happy – Neutral');
xline(0,'k--','LineWidth',1.25);
cb = colorbar; cb.Label.String = 'Power diff (a.u.)';
saveas(fh, fullfile(figureSavePath,'extra_TF_ITPC_diff_Happy_minus_Neutral_pm1s.jpg'));

% Sad - Neutral
fh = figure('Color','w','Position',[150 150 540 460]);
imagesc(tRelTF, [], diff_SN);
yticks([1:10:100])
    yticklabels(round(tfFrex(1:10:100)))
set(gca,'YDir','normal'); caxis([-.1 .1]); 
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('TF difference: Sad – Neutral');
xline(0,'k--','LineWidth',1.25);
cb = colorbar; cb.Label.String = 'Power diff (a.u.)';
saveas(fh, fullfile(figureSavePath,'extra_TF_ITPC_diff_Sad_minus_Neutral_pm1s.jpg'));


%% 4b) Gamma ROI power (±100 ms, ±5 Hz) boxplot by condition
% -------------------------------------------------------------------------
% pow_whole: events x freq x time
dt_s      = 0.1;                          % 100 ms
dt_smp    = round(dt_s * fs);
centerIdx = round((nTime + 1)/2);
tMask     = (centerIdx-dt_smp):(centerIdx+dt_smp);
tMask     = tMask(tMask>=1 & tMask<=nTime);

fMask = tfFrex >= (peakGamHz-5) & tfFrex <= (peakGamHz+5);
if ~any(fMask)
    % fallback to nearest frequency if band empty
    [~,iClosest] = min(abs(tfFrex-peakGamHz));
    fMask(iClosest) = true;
end

gammaROI = nan(nEvents,1);
for ev = 1:nEvents
    seg = pow_whole(ev, fMask, tMask);
    gammaROI(ev) = mean(seg(:),'omitnan');
end
gammaROI(burstType==0) = []; 
tmp = burstType; 
tmp(burstType == 0) = []; 
fh = figure('Color','w');
boxplot(gammaROI, tmp, 'Labels',Cnames);
ylabel('Gamma ROI power (a.u.)');
title(sprintf('Gamma power (±100 ms, %.1f±5 Hz) by condition', peakGamHz));
grid on
saveas(fh, fullfile(figureSavePath,'extra_gamma_ROI_boxplot_pm100ms_pm5Hz.jpg'));

% 
% 
% 
% 
% 
% 
% 
%     [phase, pow] = multiphasevec3(gamFrex, squeeze(macOut(chanOfInterest,:)), outDat.fs, 4,0); 
%     %what is the peak gamma frequency: 
%     figure
%     plot(gamFrex, squeeze(mean(pow, 3)), 'linewidth', 3)
%     title('channel 3 power spectrum')
%     [~, maxidx] = max(squeeze(mean(pow,3))); 
%     xline(gamFrex(maxidx))
% burstThresh = prctile(squeeze(pow(:,maxidx,:)), 90); 
%     figure
%     plot(squeeze(pow(:,maxidx,:)))
%     gamPow = squeeze(pow(:,maxidx,:)); 
%     % gamPow = (gamPow - mean(gamPow)) / std(gamPow); 
%     gamPow = smoothdata(gamPow, 'gaussian', 100); 
%     plot(gamPow)
%     yline(burstThresh)
% 
% 
%    
% 
% 
%     %get gamma burst onsets: 
%     fs = outDat.fs; 
%     gamBursts = struct; 
% 
%     gamBursts.onset = find(gamPow(1:length(gamPow)-1) < burstThresh & ...
%                     gamPow(2:length(gamPow)) > burstThresh);
% 
%     %adjust onsets to nearest oscillatory peak: 
%     gamBursts.onset = arrayfun(@(x) x - 25 + find(abs(...
%                     atan2(sin(phase(1,maxidx,x-25:x+25) - (pi/2)), ...
%                     cos(phase(1,maxidx,x-25:x+25) - (pi/2)))) == ...
%                 min(abs(atan2(sin(phase(1,maxidx,x-25:x+25) - (pi/2)), ...
%                     cos(phase(1,maxidx,x-25:x+25) - (pi/2))))) ), ...
%                     gamBursts.onset);
% 
% 
% 
%     %trim off the last couple bursts to avoid hitting the end of the recording
%     gamBursts.onset(length(gamBursts.onset)-2:length(gamBursts.onset)) = [];
%     gamBursts.onset(1:2) = [];
%     idx = find(cellfun(@(x) strcmp(x, 'badTS'), outDat.labels)); 
%     if ~isempty(idx)
%         gamBursts.bad = outDat.data(idx,gamBursts.onset);
%     else
%         gamBursts.bad = zeros(length(gamBursts.onset),1); 
%     end
%     gamBursts.onset(gamBursts.bad == 1) = []; 
%     gamBursts.bad(gamBursts.bad == 1) = []; 
%     %gamma burst lengths: 
%     gamBursts.lengths = arrayfun(@(x) find(gamPow(x:length(gamPow)-1) >burstThresh & ...
%                        gamPow(x+1:length(gamPow)) <burstThresh, 1), ...
%                        gamBursts.onset) ./ fs; 
% 
%     cycle2 = 2*(1/gamFrex(maxidx)); 
%     gamBursts.onset(gamBursts.lengths<cycle2) = [];
%     gamBursts.lengths(gamBursts.lengths<cycle2) = [];
%     gamBursts.onset(gamBursts.lengths>1) = [];
%     gamBursts.lengths(gamBursts.lengths>1) = [];
% 
%     %get the video type that each burst occurred in: 
%     gamBursts.type = typeidx(gamBursts.onset); 
% 
% 
%     %phase of respiration 
%     rspPhase = angle(hilbert(smoothdata(rspDat, 'gaussian', fs)));
%     gamBursts.rspPhase = rspPhase(gamBursts.onset); 
% 
% 
% 
%     figure
%     subplot 311
%     histogram(gamBursts.rspPhase(gamBursts.type == 1), 20);
%     title('Neutral')
%     subplot 312
%     histogram(gamBursts.rspPhase(gamBursts.type == 2), 20);
%     title('Happy')
%     subplot 313
%     histogram(gamBursts.rspPhase(gamBursts.type == 3), 20);
%     title('Sad')
% 
% %split raw ephys data into trials
% trialBuff = 4; 
%     trialDat = zeros(5, fs*2*trialBuff, length(gamBursts.onset));
%     trialRsp = zeros(fs*2*trialBuff, length(gamBursts.onset)); 
%     for trial = 1:length(gamBursts.onset)
%         curTim = gamBursts.onset(trial); 
%         curTrial = macOut(:, curTim-(fs*trialBuff):curTim+(fs*trialBuff)-1);
%         curRsp = rspDat(curTim-(fs*trialBuff):curTim+(fs*trialBuff)-1); 
%         trialDat(:,:,trial) = curTrial; 
%         trialRsp(:,trial) = curRsp; 
%     end
% 
%     badIDX = covMatClean_breathing(trialDat, ones(size(trialDat,3),1), ...
%                                     ones(size(trialDat,3),1) *...
%                                     2*fs*trialBuff); 
% 
%     trialDat(:,:,badIDX) = []; 
%     gamBursts.rspPhase(badIDX) = []; 
%     gamBursts.type(badIDX) = []; 
%     gamBursts.lengths(badIDX) = []; 
%     gamBursts.onset(badIDX) = []; 
% 
% 
% 
%     %plot an average gamma burst
%     figure
%     plot(squeeze(mean(trialDat(chanOfInterest, :, gamBursts.type==1), 3)), ...
%         'color', 'blue', 'linewidth', 2)
%     title('Neutral')
%     hold on 
%     plot(squeeze(mean(trialDat(chanOfInterest, :, gamBursts.type==2), 3)), ...
%         'color', 'green', 'linewidth', 2)
%     title('Happy')
%     plot(squeeze(mean(trialDat(chanOfInterest, :, gamBursts.type==3), 3)), ...
%         'color', 'red', 'linewidth', 2)
%     title('Sad')
% 
%     %get timefrequency decomposition for trialwise data: 
%     powDat = cell(5,2); 
%     %loop over channels
%     for ii = 3
%         ii
%         chanDat = squeeze(trialDat(ii,:,:));
%         highFrex = logspace(log10(2), log10(150), 100);
%         [phase, pow] = multiphasevec3(highFrex, chanDat', outDat.fs, 6, 0); 
%         pow = permute(pow, [3,1,2]); 
%         powDat{ii,1} = pow; 
%         % getChanMultiTF(...
%             % downDat, highFrex, 1000, .001:.001:2, 5);
%         pow = arrayfun(@(x) myChanZscore(pow(:,:,x)), 1:size(pow,3), 'UniformOutput',false ); %z-score
%         pow = cell2mat(pow); %organize
%         pow = reshape(pow, size(pow,1), size(pow,2)/length(highFrex), []); %organize
% 
% 
%         powDat{ii,2} = pow; 
% 
% 
%     end
% 
% 
% 
% 
%     figure
%     for ii = 1:3
%         subplot(3,1,ii)
%         itpc = abs(mean(exp(1i * phase(gamBursts.type==ii,:,:)), 1)); 
%         imagesc(squeeze(itpc))
%         caxis([0,.2])
%     end
% 
% 
%     figure('Color', 'w')  % white background
% 
%     timeAxis = 800:1200;                % original time vector
%     relTime = timeAxis - 1000;          % convert to relative time: -200 to +200
%     xTickVals = [800 1000 1200];
% 
%     % test = sum(pow, 1);
%     % badIDX = squeeze(sum(test,3)) > 0;
% 
%     tmp = powDat{3,2};
%     tmpbursts = gamBursts; 
%     tmpbursts.type = tmpbursts.type; 
% 
%     for i = 1:3
%         subplot(1,3,i)
% 
%         % Get average power for condition i
%         plotMat = squeeze(mean(tmp(:, tmpbursts.type==i, :), 2));
% 
%         imagesc(plotMat')  % use relTime on x-axis
%         set(gca, 'YDir', 'normal')
%         caxis([-5, 5])
% 
%         % Set custom ticks and labels
%         % xticks(xTickVals)
%         % xticklabels([-400, 0, 400])
%         % xlim([600 1400])
%         yticks([10:10:100])
%         yticklabels(round(highFrex([10:10:100])))
% 
%         % Dashed, thick vertical line at 0
%         xline(1000, 'k--', 'LineWidth', 2)
% 
%         % Title
%         switch i
%             case 1, title('Neutral')
%             case 2, title('Happy')
%             case 3, title('Sad')
%         end
% 
%         % Common labels
%         if i == 1
%             % ylabel('speed of electrical brain waves (lower to higher)')
%         end
%         % xlabel('Time relative to olfactory bulb gamma oscillation activity (milliseconds)')
%     end
% 




end

end

%% ============================== HELPERS ================================
function k = argmin_phase_distance(phiVec, target)
% Return index (1-based within phiVec) minimizing circular distance to target
d = angle(exp(1i*(phiVec - target))); % wrap to [-pi pi]
[~, k] = min(abs(d));
end

function shaded_sem(x, m, s, col, a)
% Plot mean +/- SEM as shaded patch
x = x(:)'; m = m(:)'; s = s(:)';
xx = [x, fliplr(x)];
yy = [m+s, fliplr(m-s)];
ph = patch(xx, yy, col, 'FaceAlpha', a, 'EdgeColor','none'); %#ok<NASGU>
end

function [ctrs, dens] = freqpoly(vals, nbins)
% Return mid-bin centers and normalized counts for a frequency polygon
if isempty(vals), ctrs = []; dens = []; return; end
[cts, edges] = histcounts(vals, nbins, 'Normalization','pdf');
ctrs = edges(1:end-1) + diff(edges)/2;
dens = cts;
end

function make_phase_polygons_polar(phaseVals, condType, Cnames, savePath, baseName, nBins, myTitle)
% Plot polar frequency polygons (line-joined binned counts) per condition
Ccols = lines(3);
edges = linspace(-pi, pi, nBins+1);
ctrs  = edges(1:end-1) + diff(edges)/2;

fh = figure('Color','w');
ax = polaraxes(fh); hold(ax,'on');
for c = 1:3
    ph = wrapToPi(phaseVals(condType==c));
    if isempty(ph), continue; end
    cts = histcounts(ph, edges, 'Normalization','pdf');
    rho = [cts, cts(1)];           % close the curve
    th  = [ctrs, ctrs(1)];
    polarplot(ax, th, rho, 'LineWidth', 2, 'Color', Ccols(c,:));
end
title(ax, myTitle);
legend(Cnames, 'Location','bestoutside');
saveas(fh, fullfile(savePath, [baseName '.jpg']));
end

function make_phase_polygons_after_cycles(phPeak, onsets, condType, Tgam, Kcycles, Cnames, savePath, baseBase, outfs)
% For each K in Kcycles, sample phase at onset + K*Tgam and polar-plot per condition
fs = 1 / mean(diff((1:numel(phPeak)))); %#ok<NASGU> % (unused; we’ll compute sample offset directly)
Ccols = lines(3);
for kk = 1:numel(Kcycles)
    dt = round(Kcycles(kk) * Tgam * (1/mean(diff(linspace(0,1,numel(phPeak)))))); %#ok<NASGU>
    % robust sample offset in samples:
    sampleOffset = round(Kcycles(kk) * Tgam * outfs);
    idx = onsets + sampleOffset;
    good = idx >= 1 & idx <= numel(phPeak);
    idx = idx(good);
    types = condType(good);
    phasesK = phPeak(idx);
    myTitle = ['Gamma phase dist after ' num2str(Kcycles(kk)) ' cycles'];
    make_phase_polygons_polar(phasesK, types, Cnames, savePath, ...
        sprintf('%s_K%gcycles', baseBase, Kcycles(kk)), 32, myTitle);
end
end

function tightfig()
% Reduce figure whitespace
drawnow; set(gcf,'Units','pixels');
outerpos = get(gca,'OuterPosition'); ti = get(gca,'TightInset');
left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
set(gca,'Position',[left bottom ax_width ax_height]);
end