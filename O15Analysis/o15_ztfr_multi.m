function out = o15_ztfr_multi(sig, rsp, fs, onsets, opt)
% O15_ZTFR_MULTI  Per-sniff-type baseline-z time-frequency maps for one O15
%   session's bestMac channel, all normalized to ONE shared baseline taken
%   before the first trial-start sniff.
%
%   out = o15_ztfr_multi(sig, rsp, fs, onsets [,opt])
%     onsets : struct with sample-index vectors (behDat.finalOnset, already
%              noise-cleaned upstream) for each sniff type:
%                 .start  .free  .confirm
%     opt    : optional struct overriding defaults (baseGuardS, baseDurCapS,
%              edgeTrimMs, cycles, nfreqs, frange, padPre/Post, dispPre/Post).
%
%   METHOD (analogous to cue_ztfr_pair, but a SINGLE shared baseline instead of
%   a per-trial trialStart baseline, and THREE lockings instead of two):
%     * Single baseline window = [firstStart - baseGuard - baseDurCap,
%       firstStart - baseGuard]; firstStart = min(onsets.start). One newtimef
%       pass over that segment gives the baseline TF power; wavelet edge frames
%       (+-edgeTrimMs) are trimmed -> baseFrames [F x nBaseFrames].
%     * Per sniff type: per-trial TF power from newtimef (Morlet [3 0.8], 100
%       log freqs 2-120 Hz; |alltfX|^2), averaged over trials -> meanP [F x T].
%     * Baseline-z (REUSES myChanZscore): per frequency f, z = myChanZscore(
%       [baseFrames(f,:)'; meanP(f,:)'], [1 nBaseFrames]). Passing the trial-mean
%       as a SINGLE column makes ntrials=1, so myChanZscore returns a PLAIN z =
%       (power - baselineMean)/baselineSD (baseline-SD units). Using the same
%       single baseline for all three types means the three maps share one scale
%       and are directly comparable despite very different trial counts (start/
%       confirm ~15, free ~120) -- a SEM-z (ntrials = trial count) would make the
%       free maps spuriously "hotter".
%
%   out: .freqs .baseline(.win,.nFrames,.durS,.trimmed)
%        .start/.free/.confirm = struct(map,times,freqs,resp,respT,nFinal) | []
%        .qc.<type> = struct(nTrials,nOOB,nFinal)

    if nargin < 5, opt = struct(); end
    def = struct('cycles',[3 0.8],'nfreqs',100,'frange',[2 120], ...
                 'padPre',1750,'padPost',3750,'dispPre',1000,'dispPost',3000, ...
                 'baseGuardS',1.0,'baseDurCapS',30,'edgeTrimMs',1500);
    fn = fieldnames(def);
    for k = 1:numel(fn), if ~isfield(opt,fn{k}) || isempty(opt.(fn{k})), opt.(fn{k}) = def.(fn{k}); end, end

    cycles = opt.cycles; nfreqs = opt.nfreqs; frange = opt.frange;
    padPre = opt.padPre; padPost = opt.padPost; dispPre = opt.dispPre; dispPost = opt.dispPost;

    T  = numel(sig);
    s0 = round(-padPre/1000*fs); s1 = round(padPost/1000*fs); nF = s1 - s0 + 1;
    respT = ((s0:s1)/fs*1000);
    nbfArgs = {'freqs',frange,'nfreqs',nfreqs,'freqscale','log', ...
               'baseline',NaN,'plotersp','off','plotitc','off','verbose','off'};

    out = struct(); out.freqs = []; out.qc = struct();
    types = {'start','free','confirm'};

    % ---- decompose each sniff type (to TF power) ----
    P = struct(); rp = struct(); idxOK = struct();
    freqs = []; times = [];
    for ti = 1:numel(types)
        tp = types{ti};
        ev = []; if isfield(onsets,tp), ev = round(onsets.(tp)(:)); end
        [ep, ok] = epochMat(sig, ev, s0, s1, T);
        [epR, ~] = epochMat(rsp, ev, s0, s1, T);
        out.qc.(tp) = struct('nTrials', numel(ev), 'nOOB', sum(~ok), 'nFinal', 0);
        idxOK.(tp) = find(ok);
        if numel(idxOK.(tp)) < 3
            P.(tp) = []; rp.(tp) = [];
            continue;
        end
        [~,~,~,tt,ff,~,~,atf] = newtimef(ep(:,idxOK.(tp)), nF, [-padPre padPost], fs, cycles, nbfArgs{:});
        P.(tp) = abs(atf).^2;            % [F x time x n]
        rp.(tp) = epR(:, idxOK.(tp));    % [nF x n] raw respiration epochs
        if isempty(freqs), freqs = ff; times = tt; end
    end
    if isempty(freqs)
        % no type had enough trials
        for ti = 1:numel(types), out.(types{ti}) = []; end
        out.baseline = []; return;
    end
    out.freqs = freqs;
    dispMask = times >= -dispPre & times <= dispPost;

    % ---- single shared baseline before the first start sniff ----
    baseFrames = [];
    bwin = [NaN NaN]; btrim = false; bdur = NaN;
    if isfield(onsets,'start') && ~isempty(onsets.start)
        firstStart = round(min(onsets.start(:)));
        baseEnd   = firstStart - round(opt.baseGuardS*fs);
        baseStart = max(1, baseEnd - round(opt.baseDurCapS*fs));
        if baseEnd - baseStart + 1 >= round(0.5*fs)     % need at least ~0.5 s
            seg = double(sig(baseStart:baseEnd)); Lb = numel(seg);
            [~,~,~,bt,bf,~,~,batf] = newtimef(seg(:), Lb, [0 (Lb-1)/fs*1000], fs, cycles, nbfArgs{:});
            Pbase = abs(batf).^2;                       % [F x nbt]
            % match to the signal freq grid (identical args -> identical, but be safe)
            if numel(bf) ~= numel(freqs) || any(abs(bf(:)-freqs(:)) > 1e-6)
                Pbase = interp1(bf(:), Pbase, freqs(:), 'linear', 'extrap');
            end
            valid = bt >= opt.edgeTrimMs & bt <= (bt(end) - opt.edgeTrimMs);
            if sum(valid) < 10                          % window too short to trim: keep central 60%
                lo = round(0.2*numel(bt)); hi = round(0.8*numel(bt));
                valid = false(size(bt)); valid(max(1,lo):min(numel(bt),hi)) = true;
                btrim = false;
            else
                btrim = true;
            end
            baseFrames = Pbase(:, valid);               % [F x nBaseFrames]
            bwin = [baseStart baseEnd]; bdur = Lb/fs;
        end
    end
    out.baseline = struct('win', bwin, 'nFrames', size(baseFrames,2), ...
                          'durS', bdur, 'trimmed', btrim);
    if isempty(baseFrames)
        for ti = 1:numel(types), out.(types{ti}) = []; end
        return;
    end

    % ---- baseline-z map per sniff type (shared baseline) ----
    for ti = 1:numel(types)
        tp = types{ti};
        if isempty(P.(tp)), out.(tp) = []; continue; end
        meanP = mean(P.(tp), 3, 'omitnan');             % [F x time]
        map = zmapShared(meanP, baseFrames);            % plain z vs shared baseline
        nUsed = numel(idxOK.(tp));
        out.qc.(tp).nFinal = nUsed;
        out.(tp) = struct('map', map(:,dispMask), 'times', times(dispMask), 'freqs', freqs, ...
            'resp', meanResp(rp.(tp)), 'respT', respT, 'nFinal', nUsed);
    end
end

% ======================================================================
function map = zmapShared(meanP, baseFrames)
% Plain baseline-z of a trial-mean TF map against a single shared baseline.
% Per frequency, pass the trial-mean as ONE column to myChanZscore so ntrials=1
% and the returned z = (power - baselineMean)/baselineSD (baseline-SD units).
    [F, Tn] = size(meanP);
    nb = size(baseFrames, 2);
    map = zeros(F, Tn);
    for f = 1:F
        M = [baseFrames(f,:)'; meanP(f,:)'];            % [(nb+Tn) x 1]
        z = myChanZscore(M, [1, nb]);
        map(f,:) = z(nb+1:end)';
    end
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
