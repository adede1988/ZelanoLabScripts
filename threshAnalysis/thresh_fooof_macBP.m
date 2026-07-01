function R = thresh_fooof_macBP(od, gammaBand)
% THRESH_FOOOF_MACBP  FOOOF the macBP channels of a thresh outDat (EEGLAB -> FieldTrip).
%
%   R = thresh_fooof_macBP(od, gammaBand)
%
%   Pipes the macBP channels through EEGLAB (eeg struct -> eeglab2fieldtrip),
%   segments the continuous data, computes a Welch-style PSD via FieldTrip
%   (ft_freqanalysis mtmfft), and fits FOOOF (cfg.output fooof_aperiodic /
%   fooof_peaks; brainstorm process_fooof under the hood). Picks bestMac AFTER a
%   noise screen that drops channels rejecting > noiseThresh (0.30) of trials by the
%   relative sharp-deflection rule (thresh_noise_trials), so a noisy channel whose
%   broadband noise mimics a gamma peak can't win. Among the surviving (clean)
%   channels: highest periodic gamma peak; else highest flattened power; if ALL
%   channels are too noisy, the least-noisy one (selectionMethod 'allNoisy_leastBad').
%
%   NOTE: the only thresh-specific change from cue is the noise-screen anchor,
%   od.TTL.start (thresh's pre-sniff marker), in place of od.TTL.trialStart.
%
%   Fields: .labels .chanIdx .freq .psd .aperiodic .peaksSpec .flattened
%     .gammaDetected .peakGammaFreq .peakGammaPower .flatGammaMax
%     .apExponent .apOffset .r2 .rejRate .noiseThresh
%     .bestIdx .bestMac .bestChanIdx .selectionMethod

    if nargin < 2 || isempty(gammaBand), gammaBand = [30 58]; end

    isMac  = cellfun(@(x) contains(char(string(x)),'macBP'), od.labels);
    macIdx = find(isMac);
    macLabs = cellfun(@(x) char(string(x)), od.labels(macIdx), 'uni', 0);

    R = struct();
    R.labels = macLabs; R.chanIdx = macIdx(:)'; R.gammaBand = gammaBand;
    R.empty = isempty(macIdx);
    if R.empty, return; end

    data = double(od.data(macIdx, :));
    if any(~isfinite(data(:)))
        data = fillmissing(data, 'linear', 2);
        data = fillmissing(data, 'nearest', 2);
    end
    fs = od.fs;
    nMac = numel(macIdx);

    % ---- EEGLAB EEG struct (continuous) -> FieldTrip ----
    EEG = eeg_emptyset();
    EEG.srate  = fs;
    EEG.nbchan = nMac;
    EEG.pnts   = size(data,2);
    EEG.trials = 1;
    EEG.xmin   = 0;
    EEG.xmax   = (EEG.pnts-1)/fs;
    EEG.data   = single(data);
    for c = 1:nMac, EEG.chanlocs(c).labels = macLabs{c}; end
    EEG = eeg_checkset(EEG);
    ft = eeglab2fieldtrip(EEG, 'preprocessing', 'none');

    % ---- segment + PSD (FieldTrip) ----
    rcfg = []; rcfg.length = 2; rcfg.overlap = 0.5;
    seg = ft_redefinetrial(rcfg, ft);

    base = []; base.method='mtmfft'; base.taper='hanning';
    base.foilim=[2 58]; base.pad='nextpow2'; base.keeptrials='no';   % fit 2-58 Hz (below line noise)

    cfgP = base; cfgP.output = 'pow';
    Fpow = ft_freqanalysis(cfgP, seg);

    foof = struct('aperiodic_mode','knee','max_peaks',6, ...        % knee: capture the low-freq 1/f bend
                  'peak_width_limits',[1 12],'peak_threshold',2, ...
                  'min_peak_height',0,'power_line','inf');           % 2-58 Hz excludes 60 -> no line notch
    cfgA = base; cfgA.output = 'fooof_aperiodic'; cfgA.fooof = foof;
    Fap = ft_freqanalysis(cfgA, seg);
    cfgK = base; cfgK.output = 'fooof_peaks'; cfgK.fooof = foof;
    Fpk = ft_freqanalysis(cfgK, seg);

    R.freq      = Fpow.freq(:)';
    R.psd       = reshapeChan(Fpow.powspctrm, nMac);
    R.aperiodic = reshapeChan(Fap.powspctrm,  nMac);
    R.peaksSpec = reshapeChan(Fpk.powspctrm,  nMac);
    R.flattened = R.psd ./ R.aperiodic;            % linear ratio (>1 = above 1/f)
    fp = Fap.fooofparams;                          % per-channel struct array

    inB = R.freq >= gammaBand(1) & R.freq <= gammaBand(2);
    R.gammaDetected = false(nMac,1);
    R.peakGammaFreq = nan(nMac,1);
    R.peakGammaPower= nan(nMac,1);
    R.flatGammaMax  = nan(nMac,1);
    R.apExponent    = nan(nMac,1);
    R.apOffset      = nan(nMac,1);
    R.r2            = nan(nMac,1);

    for m = 1:nMac
        R.flatGammaMax(m) = max(10*log10(R.flattened(m, inB)));   % dB over aperiodic
        ap = fp(m).aperiodic_params;
        if ~isempty(ap), R.apOffset(m) = ap(1); R.apExponent(m) = ap(end); end
        if isfield(fp,'r_squared') && ~isempty(fp(m).r_squared), R.r2(m) = fp(m).r_squared; end
        pp = fp(m).peak_params;            % [center height bw] per peak
        if ~isempty(pp)
            ing = pp(:,1) >= gammaBand(1) & pp(:,1) <= gammaBand(2);
            if any(ing)
                R.gammaDetected(m) = true;
                sub = pp(ing,:);
                [~, mx] = max(sub(:,2));
                R.peakGammaFreq(m)  = sub(mx,1);
                R.peakGammaPower(m) = sub(mx,2);
            end
        end
    end

    % --- noise screen: per-channel trial-rejection rate (relative sharp-deflection
    %     rule), so bestMac avoids channels whose "gamma" is actually broadband noise.
    %     Anchor = od.TTL.start (thresh pre-sniff marker). ---
    NOISE_THRESH = 0.30;                          % drop channels rejecting > this fraction
    R.noiseThresh = NOISE_THRESH;
    R.rejRate = nan(nMac, 1);
    ev = [];
    try, ev = od.TTL.start; catch, end            % thresh finals; skip screen if absent
    if ~isempty(ev)
        for m = 1:nMac
            NTm = thresh_noise_trials(data(m, :), fs, ev);
            if any(NTm.ok), R.rejRate(m) = mean(NTm.noisy(NTm.ok)); end
        end
    end
    elig = ~(R.rejRate > NOISE_THRESH);           % rejRate<=thresh OR NaN(no events) => eligible
    gd   = R.gammaDetected(:);

    if any(gd & elig)                             % best periodic gamma among clean channels
        cand = R.peakGammaPower; cand(~(gd & elig)) = -inf;
        [~, bi] = max(cand); R.selectionMethod = 'periodicPeak';
    elseif any(elig)                              % no clean gamma peak -> flattened among clean
        cand = R.flatGammaMax; cand(~elig) = -inf;
        [~, bi] = max(cand); R.selectionMethod = 'flattenedFallback';
    else                                          % every channel too noisy -> least-bad
        [~, bi] = min(R.rejRate); R.selectionMethod = 'allNoisy_leastBad';
    end
    R.bestIdx = bi;
    R.bestMac = macLabs{bi};
    R.bestChanIdx = macIdx(bi);
end

function M = reshapeChan(P, nMac)
% ensure [nMac x nFreq]
    if size(P,1) ~= nMac && size(P,2) == nMac, P = P'; end
    M = P;
end
