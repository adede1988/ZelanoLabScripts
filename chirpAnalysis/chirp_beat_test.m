function sub = chirp_beat_test(E, R, C)
% CHIRP_BEAT_TEST  Envelope beating (spec 6.6) -- confirms DUAL; strong over long overlaps.
%   sub = chirp_beat_test(E, R, C)
%   Beat env = abs(hilbert(broadband)) (AMPLITUDE, F4) over the COEXISTENCE window; remove the
%   onset/offset shape with a low-order polynomial (NEVER a 10 Hz HP, rule d); spectrum of the
%   residual; modSNR at fBeat=|f_hi-f_lo|. Surrogate-gated; baseline + matched-power contrasts.
%   A clean beat (>=2 cycles, modSNR>surrogate & > both controls) -> strong-for-DUAL; absent -> inconclusive.

    fs = E.fs; tMs = E.tMs;
    fBeat = abs(R.anchors.f_hi - R.anchors.f_lo);
    bb = chirp_bbfilt(E.dataPad, fs, C.bbBand, C);
    env = abs(hilbert(bb')');                  % amplitude envelope [nTrial x nPad]
    env = env(:, E.coreIdx);
    nTrial = size(E.dataPad,1);

    trial = struct('fBeat',{},'modSNR',{},'modIndex',{},'nBeatCycles',{},'deepNullPeriodicity',{}, ...
        'modSNR_baseline',{},'modSNR_matchedPower',{},'sigSurr',{},'included',{});
    logSNR = [];
    blWin = [-700 -100];                       % baseline window (ms)

    for i = 1:nTrial
        T = struct('fBeat',fBeat,'modSNR',NaN,'modIndex',NaN,'nBeatCycles',NaN, ...
            'deepNullPeriodicity',NaN,'modSNR_baseline',NaN,'modSNR_matchedPower',NaN, ...
            'sigSurr',false,'included',false);
        cw = R.trial(i).coexistWin;
        if ~E.valid(i) || isempty(cw) || ~(fBeat >= 2), trial(i)=T; continue; end
        cwDur = (cw(2)-cw(1))/1000;
        T.nBeatCycles = cwDur*fBeat;
        if T.nBeatCycles < C.beat.minCycles, trial(i)=T; continue; end

        segMask = tMs>=cw(1) & tMs<=cw(2);
        [T.modSNR, T.modIndex, T.deepNullPeriodicity] = beatMetrics(env(i,segMask), fs, fBeat, C);
        % baseline contrast (same duration if possible)
        blMask = tMs>=blWin(1) & tMs<=blWin(2);
        if any(blMask), T.modSNR_baseline = beatMetrics(env(i,blMask), fs, fBeat, C); end
        % matched-power contrast: late third of the burst (post-transition, ~stationary)
        mp = matchedWin(tMs, cw);
        if any(mp), T.modSNR_matchedPower = beatMetrics(env(i,mp), fs, fBeat, C); end
        % surrogate significance of modSNR
        seg = bb(i, :); seg = seg(E.coreIdx); seg = seg(segMask);
        Ssur = chirp_surrogates(seg, C.surr.n, C.rngSeed);
        sm = zeros(1,C.surr.n);
        for s = 1:C.surr.n, sm(s) = beatMetrics(abs(hilbert(Ssur(s,:).')).', fs, fBeat, C); end
        T.sigSurr = T.modSNR > prctile(sm, C.surr.pctile);
        T.included = true;
        trial(i) = T;
        if isfinite(T.modSNR) && T.modSNR>0, logSNR(end+1) = log(T.modSNR); end %#ok<AGROW>
    end

    used = arrayfun(@(s) s.included, trial);
    fracSig = NaN;
    if any(used)
        sg = arrayfun(@(s) isfield(s,'sigSurr') && ~isempty(s.sigSurr) && s.sigSurr, trial(used));
        fracSig = mean(sg);
    end
    sub.params = struct('coexistRef','per-trial coexistWin','envFitOrder',C.beat.envPolyOrder, ...
        'hpCutoff',C.beat.hpCutoffHz,'surrogateN',C.surr.n,'fBeat',fBeat);
    sub.trial = trial;
    sub.summary = struct('nUsed',sum(used),'medianLogSNR',medianSafe(logSNR),'fracSig',fracSig);
end

% ================= helpers =================
function [modSNR, modIndex, nullReg] = beatMetrics(envSeg, fs, fBeat, C)
    modSNR=NaN; modIndex=NaN; nullReg=NaN;
    envSeg = envSeg(:)'; n = numel(envSeg);
    if n < round(2*fs/fBeat), return; end
    % remove onset/offset SHAPE: low-order polynomial detrend (rule d; NO 10 Hz HP)
    tt = (0:n-1)/fs;
    p = polyfit(tt, envSeg, C.beat.envPolyOrder);
    resid = envSeg - polyval(p, tt);
    if ~isempty(C.beat.hpCutoffHz)               % optional, must be < min Df
        resid = chirp_bbfilt(resid, fs, [C.beat.hpCutoffHz fs/2-1], C);
    end
    % spectrum (zero-padded for bin resolution)
    nfft = 2^nextpow2(8*n);
    Y = abs(fft(resid, nfft)).^2; fr = (0:nfft-1)/nfft*fs;
    half = 1:floor(nfft/2);
    fr = fr(half); Y = Y(half);
    [~, kb] = min(abs(fr - fBeat));
    flank = (fr >= fBeat-C.beat.flankHz & fr <= fBeat+C.beat.flankHz);
    flank(max(1,kb-1):min(numel(fr),kb+1)) = false;   % exclude the +-1 bin around fBeat
    if any(flank), modSNR = Y(kb) / mean(Y(flank)); end
    modIndex = (max(resid)-min(resid)) / max(mean(abs(envSeg)),eps);
    % deep-null periodicity: regularity of envelope-minimum spacing vs 1/fBeat
    [~, locs] = findpeaks(-envSeg, 'MinPeakDistance', max(1,round(0.5*fs/fBeat)));
    if numel(locs) >= 3
        iv = diff(locs)/fs; nullReg = 1 - std(iv)/max(mean(iv),eps);
    end
end

function mp = matchedWin(tMs, cw)
    a = cw(1) + (2/3)*(cw(2)-cw(1));
    mp = tMs >= a & tMs <= cw(2);
end

function v = medianSafe(x), if isempty(x), v=NaN; else, v=median(x); end, end
