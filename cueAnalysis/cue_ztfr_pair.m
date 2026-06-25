function out = cue_ztfr_pair(sig, rsp, fs, ts, fo)
% CUE_ZTFR_PAIR  Per-trial bootstrap-z time-frequency for both lockings of a
%   cue session's bestMac channel, sharing the trialStart baseline.
%
%   out = cue_ztfr_pair(sig, rsp, fs, ts, fo)
%     ts, fo : per-trial trialStart / finalOnset sample indices (paired).
%   NOTE (analysis3): noisy trials are removed UPSTREAM (sharp-deflection rule,
%   cue_noise_trials); the ts/fo passed here are already clean. This function no
%   longer does z-score QC — it only drops out-of-bounds epochs and z-scores.
%
%   Method: per-trial TF power from newtimef (alltfX -> |.|^2), log freqs
%   2-120 Hz, Morlet [3 0.8]; bootstrap z via myChanZscore per frequency vs the
%   -700..-200ms trialStart baseline; the finalOnset epoch is z-scored against
%   the SAME baseline by appending its per-trial baseline frames ("tack-on").
%
%   out: .freqs  .trialStart{.map,.times,.resp,.respT,.nFinal}  .finalOnset{...}  .qc

    cycles = [3 0.8]; nfreqs = 100; frange = [2 120];
    padPre = 1750; padPost = 3750; dispPre = 1000; dispPost = 3000; baseMs = [-700 -200];

    T  = numel(sig);
    s0 = round(-padPre/1000*fs); s1 = round(padPost/1000*fs); nF = s1 - s0 + 1;
    N  = numel(ts);
    respT = ((s0:s1)/fs*1000);

    [epTS, okTS] = epochMat(sig, ts, s0, s1, T);
    [epFO, okFO] = epochMat(sig, fo, s0, s1, T);
    [rpTS, ~]    = epochMat(rsp, ts, s0, s1, T);
    [rpFO, ~]    = epochMat(rsp, fo, s0, s1, T);

    nbfArgs = {'freqs',frange,'nfreqs',nfreqs,'freqscale','log', ...
               'baseline',NaN,'plotersp','off','plotitc','off','verbose','off'};

    out = struct(); out.freqs = []; out.qc = struct();
    out.qc.nTrials = N;
    out.qc.nOOB_trialStart = sum(~okTS);
    out.qc.nOOB_finalOnset = sum(~okFO);

    % ---- trialStart decomposition (also supplies the shared baseline) ----
    idxTS = find(okTS);
    if numel(idxTS) < 3, out.trialStart = []; out.finalOnset = []; return; end
    [~,~,~,times,freqs,~,~,atfTS] = newtimef(epTS(:,idxTS), nF, [-padPre padPost], fs, cycles, nbfArgs{:});
    P_TS = abs(atfTS).^2;
    out.freqs = freqs;

    baseMask = times >= baseMs(1) & times <= baseMs(2);
    dispMask = times >= -dispPre & times <= dispPost;
    nb = sum(baseMask);
    baseTS = P_TS(:, baseMask, :);

    % ---- trialStart pipeline (z-score, no QC) ----
    [mapTS, keepTS] = zpipe(P_TS, baseTS, nb);
    out.trialStart = struct('map', mapTS(:,dispMask), 'times', times(dispMask), 'freqs', freqs, ...
        'resp', meanResp(rpTS(:,idxTS(keepTS))), 'respT', respT, 'nFinal', numel(keepTS));
    out.qc.nFinal_trialStart = numel(keepTS);

    % ---- finalOnset pipeline (tack-on trialStart baseline) ----
    idxFO = find(okFO);
    [hasBase, posInTS] = ismember(idxFO, idxTS);
    idxFOuse = idxFO(hasBase);
    posInTS  = posInTS(hasBase);
    out.qc.nDrop_noBaseline_finalOnset = sum(~hasBase);

    if numel(idxFOuse) >= 3
        [~,~,~,timesF,~,~,~,atfFO] = newtimef(epFO(:,idxFOuse), nF, [-padPre padPost], fs, cycles, nbfArgs{:});
        P_FO = abs(atfFO).^2;
        baseFO = baseTS(:, :, posInTS);
        [mapFO, keepFO] = zpipe(P_FO, baseFO, nb);
        out.finalOnset = struct('map', mapFO(:,dispMask), 'times', timesF(dispMask), 'freqs', freqs, ...
            'resp', meanResp(rpFO(:,idxFOuse(keepFO))), 'respT', respT, 'nFinal', numel(keepFO));
        out.qc.nFinal_finalOnset = numel(keepFO);
    else
        out.finalOnset = [];
        out.qc.nFinal_finalOnset = 0;
    end
end

% ======================================================================
function [map, keepIdx] = zpipe(Pfull, baseFrames, nb)
% Pfull [F x T x n]; baseFrames [F x nb x n] (paired). Bootstrap-z each trial per
% frequency vs the appended baseline, then average all trials (no QC; noise was
% removed upstream). The trial-average of the per-trial bootstrap z equals the
% z-of-the-mean (mean/std(dist) are per-freq scalars).
    [F, Tn, n] = size(Pfull);
    zP = zeros(F, Tn, n);
    for f = 1:F
        sigf = reshape(Pfull(f,:,:), [Tn, n]);
        basf = reshape(baseFrames(f,:,:), [nb, n]);
        M = [basf; sigf];
        z = myChanZscore(M, [1, nb]);
        zP(f,:,:) = z(nb+1:end, :);
    end
    map = mean(zP, 3, 'omitnan');
    keepIdx = 1:n;
end

% ----------------------------------------------------------------------
function [ep, ok] = epochMat(x, ev, s0, s1, T)
    nF = s1 - s0 + 1; N = numel(ev);
    ep = nan(nF, N); ok = false(N,1);
    for i = 1:N
        e = ev(i);
        if ~isfinite(e) || e <= 0, continue; end
        a = e + s0; b = e + s1;
        if a < 1 || b > T, continue; end
        seg = x(a:b);
        if any(~isfinite(seg)), continue; end
        ep(:,i) = seg(:); ok(i) = true;
    end
end

function r = meanResp(rmat)
    if isempty(rmat), r = []; return; end
    r = mean(rmat, 2, 'omitnan')';
end
