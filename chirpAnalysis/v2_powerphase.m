function pp = v2_powerphase(T, ridgeInfo, C)
% V2_POWERPHASE  Power/phase continuity across time & frequency (spec 4, powerPhaseContinuity).
%   pp = v2_powerphase(T, ridgeInfo, C)
%
%   Per good trial, using the primary ridge:
%    (1) gammaPeak{Time,Frequency,Power}: smoothed ridge-power max in C.pp.peakWinMs.
%    (2) peakBurst{Onset,Offset,Length}, hasPeak: z-power at gammaPeakFrequency crossing z=2.
%    (3) peakPhase{On/Off}{Narrow/Wide}: divergence between the narrowband phase (ridge freq at
%        peak) and a CONSTANT-frequency reference sinusoid phase-locked at the peak.
%    (4) ridgePhase{On/Off}{Narrow/Wide}: divergence between the OBSERVED phase (narrowband at
%        the moving ridge freq) and an ESTIMATED phase accumulated from the ridge freq trajectory
%        (single-oscillator hypothesis), carried forward without recentering.
%   Full phase-difference traces saved as peakPhaseConsistency / ridgePhaseConsistency
%   [N x nCore]. All outputs NaN for noisy/invalid trials; burst fields also NaN where hasPeak==0.

    fs = T.fs; tMs = T.tMs; nCore = T.nCore; N = T.N; freqs = T.freqs(:)';
    coreIdx = T.coreIdx;
    sm = C.pp.smoothWin; thN = C.pp.threshNarrow; thW = C.pp.threshWide; k5 = C.pp.phaseAvgSamp;
    dt = 1/fs;

    nan1 = nan(N,1); nanM = nan(N,nCore);
    pp = struct('gammaPeakTime',nan1,'gammaPeakFrequency',nan1,'gammaPeakPower',nan1, ...
        'peakBurstOnset',nan1,'peakBurstOffset',nan1,'peakBurstLength',nan1,'hasPeak',nan1, ...
        'ridgeBurstOnset',nan1,'ridgeBurstOffset',nan1,'ridgeBurstLength',nan1,'burstTruncated',nan1, ...
        'peakPhaseOnsetNarrow',nan1,'peakPhaseOffsetNarrow',nan1,'peakPhaseOnsetWide',nan1,'peakPhaseOffsetWide',nan1, ...
        'ridgePhaseOnsetNarrow',nan1,'ridgePhaseOffsetNarrow',nan1,'ridgePhaseOnsetWide',nan1,'ridgePhaseOffsetWide',nan1, ...
        'peakPhaseConsistency',nanM,'ridgePhaseConsistency',nanM,'tMs',tMs, ...
        'params',struct('smoothWin',sm,'burstZ',C.pp.burstZ,'narrowBWhz',C.pp.narrowBWhz, ...
                        'threshNarrow',thN,'threshWide',thW,'phaseAvgSamp',k5,'peakWinMs',C.pp.peakWinMs));

    prF = ridgeInfo.primaryRidge.f; prP = ridgeInfo.primaryRidge.p;
    winMask = tMs >= C.pp.peakWinMs(1) & tMs <= C.pp.peakWinMs(2);
    idxWin = find(winMask);

    for i = find(T.good)'
        f_r = prF(i,:); p_r = prP(i,:);
        if all(~isfinite(p_r)) || isempty(idxWin), continue; end

        % (1) gamma peak (smoothed ridge power, max in window) -----------------------------
        pS = smoothdata(p_r, 'gaussian', sm);
        if all(~isfinite(pS(idxWin))), continue; end   % peak window fully truncated (O15 guard)
        [~, mrel] = max(pS(idxWin)); kpk = idxWin(mrel);
        gFreq = f_r(kpk);
        pp.gammaPeakTime(i)      = tMs(kpk);
        pp.gammaPeakFrequency(i) = gFreq;
        pp.gammaPeakPower(i)     = p_r(kpk);            % UNsmoothed ridge power at peak

        % (2) peak burst (z-power at gammaPeakFrequency crossing z=2) ----------------------
        [~, fk] = min(abs(freqs - gFreq));
        zser = reshape(T.zTFR(fk,:,i), 1, nCore);
        zS   = smoothdata(zser, 'gaussian', sm);
        hasPk = zS(kpk) >= C.pp.burstZ;
        pp.hasPeak(i) = double(hasPk);
        if hasPk
            on = kpk;  while on>1     && zS(on-1) >= C.pp.burstZ, on = on-1;   end
            off= kpk;  while off<nCore&& zS(off+1)>= C.pp.burstZ, off= off+1;  end
            pp.peakBurstOnset(i)  = tMs(on);
            pp.peakBurstOffset(i) = tMs(off);
            pp.peakBurstLength(i) = tMs(off) - tMs(on);
        end

        % (2b) ridge burst (smoothed PRIMARY RIDGE power crossing z=2, from gammaPeakTime) --
        if pS(kpk) >= C.pp.burstZ
            on = kpk;  while on>1     && pS(on-1) >= C.pp.burstZ, on = on-1;   end
            off= kpk;  while off<nCore&& pS(off+1)>= C.pp.burstZ, off= off+1;  end
            pp.ridgeBurstOnset(i)  = tMs(on);
            pp.ridgeBurstOffset(i) = tMs(off);
            pp.ridgeBurstLength(i) = tMs(off) - tMs(on);
        end

        % (3) peak-based phase continuity vs constant-frequency reference ------------------
        nb  = chirp_bbfilt(T.padData(i,:), fs, gFreq + [-1 1]*C.pp.narrowBWhz/2, C);
        anb = hilbert(nb(:)).';                           % hilbert needs a column; return row
        phiObs = angle(anb(coreIdx));                    % 1 x nCore
        tRelS  = (tMs - tMs(kpk))/1000;                  % s from peak
        phiRef = phiObs(kpk) + 2*pi*gFreq*tRelS;         % constant-freq reference phase
        dPeak  = wrapToPi(phiObs - phiRef);
        pp.peakPhaseConsistency(i,:) = dPeak;
        adPeak = circSmoothAbs(dPeak, k5);
        [pp.peakPhaseOnsetNarrow(i), pp.peakPhaseOffsetNarrow(i)] = crossings(adPeak, kpk, thN, tMs);
        [pp.peakPhaseOnsetWide(i),   pp.peakPhaseOffsetWide(i)]   = crossings(adPeak, kpk, thW, tMs);

        % (4) ridge-based phase progression (accumulated single-oscillator estimate) -------
        fInt = max(C.ridge.band(1), floor(min(f_r))) : min(C.ridge.band(2), ceil(max(f_r)));
        if isempty(fInt), fInt = round(gFreq); end
        phiBank = zeros(numel(fInt), nCore);             % narrowband phase per integer freq
        for q = 1:numel(fInt)
            nbq = chirp_bbfilt(T.padData(i,:), fs, fInt(q) + [-1 1]*C.pp.narrowBWhz/2, C);
            aq  = hilbert(nbq(:)).';                      % hilbert needs a column; return row
            phiBank(q,:) = angle(aq(coreIdx));
        end
        rf = min(max(round(f_r), fInt(1)), fInt(end));
        if numel(fInt) >= 2
            obsIdx = interp1(fInt, 1:numel(fInt), rf, 'nearest', 'extrap');
        else
            obsIdx = ones(1, nCore);
        end
        obsIdx = min(max(round(obsIdx),1), numel(fInt));
        phiObsR = phiBank(sub2ind(size(phiBank), obsIdx, 1:nCore));   % observed phase at ridge freq
        % accumulate estimated phase outward from the peak
        phiEst = nan(1,nCore); phiEst(kpk) = phiObsR(kpk);
        for t = kpk+1:nCore,  phiEst(t) = phiEst(t-1) + 2*pi*f_r(t-1)*dt; end
        for t = kpk-1:-1:1,   phiEst(t) = phiEst(t+1) - 2*pi*f_r(t+1)*dt; end
        dRidge = wrapToPi(phiEst - phiObsR);
        pp.ridgePhaseConsistency(i,:) = dRidge;
        adR = circSmoothAbs(dRidge, k5);
        [pp.ridgePhaseOnsetNarrow(i), pp.ridgePhaseOffsetNarrow(i)] = crossings(adR, kpk, thN, tMs);
        [pp.ridgePhaseOnsetWide(i),   pp.ridgePhaseOffsetWide(i)]   = crossings(adR, kpk, thW, tMs);
    end
    if isfield(T,'burstTruncated'), pp.burstTruncated = double(T.burstTruncated); end
end

% ---- circular running-mean magnitude of a wrapped phase-difference series ----
function ad = circSmoothAbs(d, k)
    dd = angle(movmean(exp(1i*d), k));   % circular running mean (robust near +-pi)
    ad = abs(dd);
end

% ---- first threshold crossing walking back (onset) / forward (offset) from kpk ----
function [onMs, offMs] = crossings(ad, kpk, thr, tMs)
    nC = numel(ad);
    on = kpk;  while on>1  && ad(on-1)  <= thr, on = on-1;   end
    off= kpk;  while off<nC && ad(off+1) <= thr, off= off+1;  end
    onMs = tMs(on); offMs = tMs(off);
end
