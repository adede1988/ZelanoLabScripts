function out = thresh_ztfr_conds(sig, rsp, fs, startByCond, foByCond)
% THRESH_ZTFR_CONDS  Per-odor-condition bootstrap-z time-frequency maps for one
%   thresh session's bestMac channel. The core scientific change from cue: instead
%   of two lockings (trialStart / finalOnset) sharing one baseline, we have THREE
%   odor conditions (low / med / air) that are ALL finalOnset-locked, each z-scored
%   against ITS OWN trials' pre-`start` baseline (the same baseline DEFINITION as
%   cue's finalOnset branch, applied per condition).
%
%   out = thresh_ztfr_conds(sig, rsp, fs, startByCond, foByCond)
%     startByCond : struct .low .med .air, each a vector of TTL.start sample indices
%                   (the per-trial baseline anchor), clean trials only.
%     foByCond    : struct .low .med .air, each a vector of finalOnset sample indices
%                   for that condition's clean trials, PAIRED element-wise with the
%                   matching entries of startByCond (so each condition trial finds
%                   its own pre-`start` baseline frames).
%   NOTE: noisy trials are removed UPSTREAM (relative sharp-deflection rule,
%   thresh_noise_trials); the vectors passed here are already clean. This function
%   only drops out-of-bounds epochs and z-scores.
%
%   Method (identical engine to cue's finalOnset branch, run three times):
%   1. Per condition: newtimef (log freqs 2-120 Hz, Morlet [3 0.8], baseline NaN ->
%      raw power |alltfX|^2) on the condition's `start` epochs (for the baseline
%      frames) and its `finalOnset` epochs.
%   2. Baseline = -700..-200 ms of the trials' `start` epochs; the per-trial baseline
%      frames are tacked onto the front of each finalOnset trial and myChanZscore'd
%      per frequency (reuse the `zpipe` subfunction verbatim). Average trials -> map.
%   3. Respiration: mean of that condition's finalOnset-locked `rsp` epochs (overlay).
%
%   out: .freqs .low/.med/.air = struct(map,times,resp,respT,nFinal) | [] (if <3
%        clean trials) .qc.<cond> = struct(nTrials, nOOB, nFinal)

    cycles = [3 0.8]; nfreqs = 100; frange = [2 120];
    padPre = 1750; padPost = 3750; dispPre = 1000; dispPost = 3000; baseMs = [-700 -200];

    T  = numel(sig);
    s0 = round(-padPre/1000*fs); s1 = round(padPost/1000*fs); nF = s1 - s0 + 1;
    respT = ((s0:s1)/fs*1000);
    nbfArgs = {'freqs',frange,'nfreqs',nfreqs,'freqscale','log', ...
               'baseline',NaN,'plotersp','off','plotitc','off','verbose','off'};

    conds = {'low','med','air'};
    out = struct(); out.freqs = []; out.qc = struct();
    freqs = []; times = []; baseMask = []; dispMask = [];

    for ci = 1:numel(conds)
        c = conds{ci};
        ts = []; fo = [];
        if isfield(startByCond, c), ts = round(startByCond.(c)(:)); end
        if isfield(foByCond, c),    fo = round(foByCond.(c)(:)); end
        out.qc.(c) = struct('nTrials', numel(fo), 'nOOB', 0, 'nFinal', 0);
        out.(c) = [];
        if numel(fo) < 3 || numel(ts) ~= numel(fo), continue; end

        [epTS, okTS] = epochMat(sig, ts, s0, s1, T);
        [epFO, okFO] = epochMat(sig, fo, s0, s1, T);
        [rpFO, ~]    = epochMat(rsp, fo, s0, s1, T);
        out.qc.(c).nOOB = sum(~okFO);

        idxTS = find(okTS);
        if numel(idxTS) < 3, continue; end
        [~,~,~,tt,ff,~,~,atfTS] = newtimef(epTS(:,idxTS), nF, [-padPre padPost], fs, cycles, nbfArgs{:});
        P_TS = abs(atfTS).^2;
        if isempty(freqs)
            freqs = ff; times = tt; out.freqs = freqs;
            baseMask = times >= baseMs(1) & times <= baseMs(2);
            dispMask = times >= -dispPre & times <= dispPost;
        end
        nb = sum(baseMask);
        baseTS = P_TS(:, baseMask, :);

        % pair each valid finalOnset trial with its own (valid) start baseline
        idxFO = find(okFO);
        [hasBase, posInTS] = ismember(idxFO, idxTS);
        idxFOuse = idxFO(hasBase);
        posInTS  = posInTS(hasBase);
        if numel(idxFOuse) < 3, continue; end

        [~,~,~,~,~,~,~,atfFO] = newtimef(epFO(:,idxFOuse), nF, [-padPre padPost], fs, cycles, nbfArgs{:});
        P_FO = abs(atfFO).^2;
        baseFO = baseTS(:, :, posInTS);
        [mapFO, keepFO] = zpipe(P_FO, baseFO, nb);
        out.(c) = struct('map', mapFO(:,dispMask), 'times', times(dispMask), 'freqs', freqs, ...
            'resp', meanResp(rpFO(:,idxFOuse(keepFO))), 'respT', respT, 'nFinal', numel(keepFO));
        out.qc.(c).nFinal = numel(keepFO);
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
