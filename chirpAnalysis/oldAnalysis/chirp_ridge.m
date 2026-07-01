function R = chirp_ridge(tfr, E, C)
% CHIRP_RIDGE  Ridge extraction, peeling, transition/coexistence windows, dip stat, anchors.
%   R = chirp_ridge(tfr, E, C)
%
%   Per valid trial: primary tfridge ridge fhat(t)+ridgePower(t) in 25-58 Hz; per-trial
%   transitionWin (from THIS trial's own smoothed fhat extrema -- F2) + coexistWin (burst span);
%   ridge-power-dip statistic (burst-relative flanks -- F7); peeled second ridge (+-peelHz).
%
%   2nd-ridge significance (F6, cost-aware): the null distribution of the normalized peeled/primary
%   ridge-energy RATIO is a per-CHANNEL property, so it is estimated ONCE from phase-randomized
%   surrogates of a few sampled trials (NOT per trial) -> ~20x fewer FASLTs than per-trial surrogates.
%   nRidgesSig = 1 + (trial ratio2 > 95th-pctile null ratio2). All windows in MILLISECONDS rel onset.

    f   = tfr.freqs(tfr.bandMask); tMs = tfr.tMs; fs = tfr.fs; nT = numel(tMs);
    pen = C.ridge.penalty;
    peelW = max(1, round(C.ridge.peelHalfHz / median(diff(f))));
    smW  = max(1, round(C.ridge.fhatSmoothMs/1000*fs));
    nTrial = size(tfr.power, 3);

    trial = repmat(blankTrial(), nTrial, 1);
    ratio2 = nan(nTrial,1);
    for i = 1:nTrial
        if ~tfr.valid(i), continue; end
        tfm = tfr.power(tfr.bandMask, :, i);
        if all(~isfinite(tfm(:))) || all(tfm(:)==0), continue; end
        T = blankTrial();
        [fhat, iR] = tfridge(tfm, f, pen, 'NumRidges',1);
        fhat = fhat(:)'; iR = iR(:)';
        rp = tfm(sub2ind(size(tfm), iR, 1:nT));
        fhatS = movmedian(fhat, smW);

        post = tMs > 0 & tMs <= C.burst.searchMs(2);
        [bm, b0, b1] = burstSpan(rp, post);
        T.fhat = fhat; T.ridgePower = rp;
        if isempty(b0), trial(i) = T; continue; end
        T.coexistWin = [tMs(b0) tMs(b1)];
        T.included = isfinite(rp(b0));

        bi = b0:b1; third = max(1, floor(numel(bi)/3));
        fHiTr = median(fhatS(bi(1:third)), 'omitnan');
        fLoTr = median(fhatS(bi(end-third+1:end)), 'omitnan');
        if (fHiTr - fLoTr) >= 2*C.ridge.transEdgeHz
            kS = find(fhatS(bi) <= fHiTr - C.ridge.transEdgeHz, 1, 'first');
            kE = find(fhatS(bi) <= fLoTr + C.ridge.transEdgeHz, 1, 'first');
            if ~isempty(kS) && ~isempty(kE) && kE > kS
                T.transitionWin = [tMs(bi(kS)) tMs(bi(kE))];
            end
        end
        T.powerDipStat = dipStat(rp, tMs, bi, T.transitionWin, fs);

        % peeled second ridge
        tfm2 = peel(tfm, iR, peelW);
        [fhat2, iR2] = tfridge(tfm2, f, pen, 'NumRidges',1);
        fhat2 = fhat2(:)'; rp2 = tfm2(sub2ind(size(tfm2), iR2(:)', 1:nT));
        T.fhat2 = fhat2; T.ridgePower2 = rp2;

        win = T.transitionWin; if isempty(win), win = T.coexistWin; end
        wm = tMs >= win(1) & tMs <= win(2);
        e1 = sum(rp(wm),'omitnan'); e2 = sum(rp2(wm),'omitnan');
        ratio2(i) = e2 / max(e1, eps);
        T.bm = bm;
        trial(i) = T;
    end

    % ---- per-channel null ratio2 (sampled trials x surrogates), F6 cost-aware ----
    thr95 = inf;
    vi = find(tfr.valid(:)' & arrayfun(@(t) ~isempty(t.coexistWin), trial(:)'));
    if ~isempty(vi)
        nNull = min(getf(C.ridge,'nNullTrials',3), numel(vi));
        nSur  = getf(C.ridge,'nSurr',60);
        samp = vi(round(linspace(1, numel(vi), nNull)));
        nullR = [];
        for s = samp
            cw = trial(s).coexistWin; wm = tMs>=cw(1) & tMs<=cw(2);
            Ssur = chirp_surrogates(E.dataPad(s,:), nSur, C.rngSeed + s);
            for q = 1:nSur
                wt = chirp_faslt_apply(Ssur(q,:), tfr.bank);
                tm = wt(tfr.bandMask, E.coreIdx);
                [~, iRs] = tfridge(tm, f, pen, 'NumRidges',1);
                rps = tm(sub2ind(size(tm), iRs(:)', 1:nT));
                tm2 = peel(tm, iRs(:)', peelW);
                [~, iRs2] = tfridge(tm2, f, pen, 'NumRidges',1);
                rps2 = tm2(sub2ind(size(tm2), iRs2(:)', 1:nT));
                nullR(end+1) = sum(rps2(wm),'omitnan')/max(sum(rps(wm),'omitnan'),eps); %#ok<AGROW>
            end
        end
        if ~isempty(nullR), thr95 = prctile(nullR, C.surr.pctile); end
    end

    % ---- assign nRidgesSig + ridgeOverlap ----
    for i = 1:nTrial
        if ~trial(i).included, continue; end
        sig2 = ratio2(i) > thr95;
        trial(i).nRidgesSig = 1 + double(sig2);
        if sig2 && ~isempty(trial(i).bm)
            rp2 = trial(i).ridgePower2; bm = trial(i).bm;
            ov = (rp2 > 0.5*max(rp2(bm))) & bm & (abs(trial(i).fhat2 - trial(i).fhat) > C.chirplet.dualDf);
            trial(i).ridgeOverlap = sum(ov)/max(1,sum(bm));
        else
            trial(i).ridgeOverlap = 0;
        end
    end
    trial = rmfield(trial, 'bm');

    R.trial   = trial;
    R.tMs     = tMs; R.freqsB = f;
    R.ratio2  = ratio2; R.nullThr95 = thr95;
    R.anchors = chirp_ridge_anchors(tfr, C);
end

% ===================== helpers =====================
function T = blankTrial()
    T = struct('fhat',[],'ridgePower',[],'fhat2',[],'ridgePower2',[],'nRidgesSig',1, ...
        'ridgeOverlap',0,'powerDipStat',NaN,'transitionWin',[],'coexistWin',[], ...
        'burstTruncated',false,'included',false,'bm',[]);
end
function [bm, b0, b1] = burstSpan(rp, post)
    bm = false(size(rp)); b0 = []; b1 = [];
    idx = find(post); if isempty(idx), return; end
    [pk, krel] = max(rp(idx)); if ~(pk>0), return; end
    thr = 0.5*pk; kpk = idx(krel);
    a = kpk; while a>1 && post(a-1) && rp(a-1) > thr, a = a-1; end
    b = kpk; while b<numel(rp) && post(b+1) && rp(b+1) > thr, b = b+1; end
    bm(a:b) = true; b0 = a; b1 = b;
end
function d = dipStat(rp, tMs, bi, transWin, fs)
    d = NaN; if isempty(transWin), return; end
    inTr = tMs >= transWin(1) & tMs <= transWin(2); if ~any(inTr), return; end
    flankN = max(1, round(0.15*fs));
    pre = bi(1):min(bi(1)+flankN-1, bi(end));
    post = max(bi(1), bi(end)-flankN+1):bi(end);
    fmax = mean([max(rp(pre)), max(rp(post))], 'omitnan');
    if ~(fmax > 0), return; end
    d = min(rp(inTr)) / fmax;
end
function tfm = peel(tfm, iR, W)
    nT = size(tfm,2); nF = size(tfm,1);
    for t = 1:nT, tfm(max(1,iR(t)-W):min(nF,iR(t)+W), t) = 0; end
end
function v = getf(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
