function sub = chirp_temporal_decoup(E, R, C)
% CHIRP_TEMPORAL_DECOUP  Cross-trial latency decoupling (spec 6.9) -- pairs with phase; leans dual.
%   sub = chirp_temporal_decoup(E, R, C)
%   Band-limited POWER envelopes (FIR-hilbert, F4) at f_hi and f_lo on the OB channel; per-trial
%   peak latency of each + their difference. One chirping oscillator -> rigid coupling
%   (corr(latHi,latLo)->1, low var(latDiff)); two independently-timed generators -> decoupled.

    fs = E.fs; tMs = E.tMs;
    fHi = R.anchors.f_hi; fLo = R.anchors.f_lo; bh = C.temporal.bandHalfHz;
    nTrial = size(E.dataPad,1);
    % early exit (anchors collapsed/non-finite): return a BLANK trial PER trial (sized to
    % nTrial), not an empty struct -- appendTrialCSVs loops i=1:nTr and would index past an
    % empty array (crashes on sessions whose ridge anchors fall within 1 Hz).
    if ~(isfinite(fHi)&&isfinite(fLo)) || abs(fHi-fLo) < 1
        bt = struct('latHi',NaN,'latLo',NaN,'latDiff',NaN,'included',false);
        sub = struct('trial', repmat(bt,1,nTrial), ...
            'summary', struct('corrHiLo',NaN,'varLatDiff',NaN,'nUsed',0)); return;
    end
    hiP = bandPow(E.dataPad, fs, [fHi-bh fHi+bh], C); hiP = hiP(:,E.coreIdx);
    loP = bandPow(E.dataPad, fs, [fLo-bh fLo+bh], C); loP = loP(:,E.coreIdx);
    smW = max(1, round(C.temporal.smoothMs/1000*fs));
    hiP = movmean(hiP, smW, 2); loP = movmean(loP, smW, 2);

    trial = struct('latHi',{},'latLo',{},'latDiff',{},'included',{});
    LH=[]; LL=[]; LD=[];
    for i = 1:nTrial
        T = struct('latHi',NaN,'latLo',NaN,'latDiff',NaN,'included',false);
        cw = R.trial(i).coexistWin;
        if ~E.valid(i) || isempty(cw), trial(i)=T; continue; end
        m = tMs>=cw(1) & tMs<=cw(2); idx = find(m);
        if isempty(idx), trial(i)=T; continue; end
        [pkH,kH] = max(hiP(i,idx)); [pkL,kL] = max(loP(i,idx));
        % existence: both bands must have positive power in the window
        if ~(pkH>0 && pkL>0), trial(i)=T; continue; end
        T.latHi = tMs(idx(kH)); T.latLo = tMs(idx(kL)); T.latDiff = T.latLo - T.latHi;
        T.included = true; trial(i)=T;
        LH(end+1)=T.latHi; LL(end+1)=T.latLo; LD(end+1)=T.latDiff; %#ok<AGROW>
    end
    corrHiLo = NaN;
    if numel(LH) >= 3, c = corrcoef(LH,LL); corrHiLo = c(1,2); end
    sub.trial = trial;
    sub.summary = struct('corrHiLo',corrHiLo,'varLatDiff',varSafe(LD),'nUsed',numel(LH));
end

function P = bandPow(x, fs, band, C)
    band(1) = max(band(1), 1); band(2) = min(band(2), fs/2-1);
    y = chirp_bbfilt(x, fs, band, C);
    P = abs(hilbert(y')').^2;
end
function v = varSafe(x), if numel(x)<2, v=NaN; else, v=var(x); end, end
