function [plr, meanPeakTFR, relFreq, relTimeMs, centStack] = v2_peaklocked(T, pp, C)
% V2_PEAKLOCKED  Peak-aligned re-analysis (spec chirpAnalysis_3.md sec 4).
%   [plr, meanPeakTFR, relFreq, relTimeMs] = v2_peaklocked(T, pp, C)
%
%   For each good trial: re-window +-1000 ms around gammaPeakTime (re-cut from the already
%   padded onset-locked epoch, which carries +-1.5 s pad -> a peak in [0,700]ms leaves room for
%   the +-1000 ms window plus FASLT pad). Re-run FASLT (WIDER range 15-75 for +-10 Hz headroom
%   around the peak freq) + baseline z-score (-1000..-500 ms rel peak) + primary ridge (25-58).
%
%   Returns:
%     plr          : peakLockedRidge struct (.f .p [N x nPk], .tRelMs .freqs .band .note)
%     meanPeakTFR  : subject mean peak-aligned TFR [nRelF x nPk], freq axis CENTERED on each
%                    trial's gammaPeakFrequency (+-10 Hz), time centered on gammaPeakTime (+-1000ms)
%     relFreq [1 x nRelF] (Hz rel peak), relTimeMs [1 x nPk] (ms rel peak)

    fs = C.fs; N = T.N; tPadMs = T.tPadMs;
    coreHalf = round(1.0*fs);           % +-1000 ms core
    padHalf  = round(C.padSec*fs);      % FASLT pad
    nPk = 2*coreHalf + 1;
    relTimeMs = (-coreHalf:coreHalf)/fs*1000;

    % wider FASLT bank for +-10 Hz headroom around the peak frequency
    Cw = C; Cw.faslt.range = [15 75]; Cw.faslt.Nf = 61;
    bankW = chirp_faslt_bank(fs, Cw);
    fW = bankW.F(:)';                    % 15..75 Hz
    bandIdx = find(fW >= C.ridge.band(1) & fW <= C.ridge.band(2));
    fBand = fW(bandIdx);
    relFreq = -10:1:10; nRelF = numel(relFreq);

    % baseline indices (-1000..-500 ms rel peak) on the core grid
    bIdx = find(relTimeMs >= -1000 & relTimeMs <= -500); basePeriod = [bIdx(1) bIdx(end)];

    plF = nan(N,nPk); plP = nan(N,nPk);
    zStack = nan(61, nPk, N);           % per-trial peak-aligned z-TFR (absolute freq 15-75)
    centStack = nan(nRelF, nPk, N);     % per-trial freq-centered z-TFR (rel freq -10..10)
    gi = [];

    for i = find(T.good)'
        gpt = pp.gammaPeakTime(i); gpf = pp.gammaPeakFrequency(i);
        if ~isfinite(gpt) || ~isfinite(gpf), continue; end
        kp = round((gpt - tPadMs(1))/1000*fs) + 1;      % peak index in padded epoch
        a = kp - coreHalf - padHalf; b = kp + coreHalf + padHalf;
        if a < 1 || b > numel(tPadMs), continue; end
        seg = T.padData(i, a:b);
        if any(~isfinite(seg)), continue; end
        wt = chirp_faslt_apply(seg, bankW);              % [61 x (nPk+2*padHalf)]
        core = (padHalf+1):(padHalf+nPk);
        zStack(:,:,i) = wt(:, core);                     % trim pad -> [61 x nPk]  (POWER; z below)
        gi(end+1) = i; %#ok<AGROW>
    end

    % per-frequency baseline z-score across good trials (myChanZscore), then ridge + centering
    if ~isempty(gi)
        for f = 1:61
            dataF = reshape(zStack(f,:,gi), nPk, numel(gi));
            zStack(f,:,gi) = reshape(myChanZscore(dataF, basePeriod), 1, nPk, numel(gi));
        end
        for i = gi
            Z = zStack(:,:,i);
            % primary ridge on the peak-aligned z-TFR
            Zb = Z(bandIdx,:);
            [fr, ir] = tfridge(max(Zb,0), fBand, C.ridge.penalty, 'NumRidges',1);
            plF(i,:) = fr(:)';
            plP(i,:) = Zb(sub2ind(size(Zb), ir(:)', 1:nPk));
            % freq-centered slice: interp rows at gammaPeakFrequency + relFreq
            absF = pp.gammaPeakFrequency(i) + relFreq;
            centStack(:,:,i) = interp1(fW, Z, absF, 'linear', NaN);   % [nRelF x nPk]
        end
    end

    meanPeakTFR = mean(centStack(:,:,gi), 3, 'omitnan');
    plr = struct('f',plF,'p',plP,'tRelMs',relTimeMs,'freqs',fW,'band',C.ridge.band, ...
        'note','peak-aligned (+-1000ms around gammaPeakTime) FASLT 15-75, z-score, tfridge primary');
end
