function out = chirp_chirplet(segCore, tMs, transWin, C)
% CHIRP_CHIRPLET  MPACT adaptive-chirplet trajectory-GEOMETRY discriminator (spec 6.7, rule c).
%   out = chirp_chirplet(segCore, tMs, transWin, C)
%   segCore : 1 x nCore RAW (unfiltered) core-trimmed epoch of the OB channel
%   tMs     : 1 x nCore time axis (ms rel finalOnset) matching segCore
%   transWin: [t0 t1] ms (per-trial transition window from chirp_ridge) -- fit ONLY here
%   C       : chirp_config
%
%   Bandpass 25-58, crop to the transition window, feed the ANALYTIC signal to
%   mp_adapt_chirplets (MPEM, verbose='no', NOT mp_act_signal/mle_adapt_chirplets). Classify by
%   trajectory GEOMETRY: single = atoms chain into one descending trajectory; dual = >=2
%   near-constant-freq atoms at distinct fc that temporally COEXIST. dBIC is corroborative only.

    fs = C.fs;
    out = blank();
    if isempty(transWin) || any(isnan(transWin)), out.reason='noTransWin'; return; end

    bb = chirp_bbfilt(segCore, fs, C.bbBand, C);
    % cost guard: cap the fit window length (the chirp action is early; MPEM cost grows fast with N)
    w2 = transWin(2);
    if isfield(C.chirplet,'maxWinMs') && ~isempty(C.chirplet.maxWinMs)
        w2 = min(transWin(2), transWin(1) + C.chirplet.maxWinMs);
    end
    idx = tMs >= transWin(1) & tMs <= w2;
    crop = bb(idx); crop = detrend(crop(:)' - mean(crop));
    Ntr = numel(crop);
    if Ntr < C.chirplet.minTransSamp, out.reason='shortTransition'; return; end

    xa = hilbert(crop); xa = xa(:);
    Q = C.chirplet.Q;
    [P, res] = mpfit(xa, Q, C);
    if isempty(P), out.reason='fitFailed'; return; end

    % convert O'Neill -> physical units
    A = P(:,1); tc = P(:,2); fc = P(:,3); cr = P(:,4); d = P(:,5);
    tc_s   = (tc-1)/fs + transWin(1)/1000;
    fc_Hz  = fc*fs/(2*pi);
    chirpHzs = cr*fs^2/(2*pi);
    sig_s  = d/fs;
    amp    = abs(A);
    ph     = angle(A);
    Etot   = sum(abs(xa).^2);
    eFrac  = (amp.^2)/max(Etot,eps);

    % atom significance: 'fixed' energy-fraction threshold (fast batch default; Phase-0 calibrated)
    % or per-trial surrogate single-atom energy fraction ('surrogate', slow -- Phase-0 calibration).
    if isfield(C.chirplet,'sigMode') && strcmp(C.chirplet.sigMode,'surrogate')
        if isfield(C.chirplet,'nSurr'), nS = C.chirplet.nSurr; else, nS = C.surr.n; end
        Ssur = chirp_surrogates(crop, nS, C.rngSeed);
        surrFrac = zeros(1,nS);
        for s = 1:nS
            xs = hilbert(Ssur(s,:)); xs = xs(:);
            Ps = mpfit(xs, 1, C);
            if ~isempty(Ps), surrFrac(s) = abs(Ps(1,1))^2/max(sum(abs(xs).^2),eps); end
        end
        thr = prctile(surrFrac, C.surr.pctile);
    else
        thr = C.chirplet.atomEnergyThr;
    end
    out.atomThr = thr;
    keep = (eFrac > thr) & (fc_Hz >= C.chirplet.atomKeepBand(1)) & (fc_Hz <= C.chirplet.atomKeepBand(2));

    atoms = [tc_s(:) fc_Hz(:) chirpHzs(:) sig_s(:) amp(:) ph(:) eFrac(:)];
    atoms = atoms(keep, :);
    [~, o] = sort(atoms(:,1)); atoms = atoms(o,:);   % by tc
    nSig = size(atoms,1);
    out.atoms = atoms; out.nAtomsSig = nSig;

    % trajectory geometry
    [traj, nChain, nPar] = geometry(atoms, C);
    out.trajConnectivity = traj;
    out.classification = classify(atoms, nChain, nPar, C);

    % chirp-reality (Model-A single atom free c) + corroborative dBIC from res
    out.c1 = NaN; out.errA = NaN; out.errB = NaN; out.dBIC = NaN;
    if numel(res) >= 3
        errA = res(2)^2; errB = res(3)^2;       % residual energy after 1 / 2 atoms (unconstrained)
        meanFc = mean(fc_Hz(isfinite(fc_Hz)));
        Neff = max(2, meanFc*Ntr/fs);           % #gamma cycles (NOT raw samples)
        bicA = Neff*log(max(errA,eps)/Neff) + 5*log(Neff);
        bicB = Neff*log(max(errB,eps)/Neff) + 10*log(Neff);
        out.errA = errA; out.errB = errB; out.dBIC = bicA - bicB;
    end
    % c1 = chirp rate of the dominant (highest-energy) kept atom
    if nSig >= 1, [~,im] = max(atoms(:,5)); out.c1 = atoms(im,3); end

    % f1/f2/tSep/overlap descriptors
    [out.f1,out.f2,out.tSep,out.overlap] = descriptors(atoms, out.classification);
    out.included = true; out.reason = '';
end

% ================= helpers =================
function out = blank()
    out = struct('atoms',[],'nAtomsSig',0,'trajConnectivity',NaN,'classification','ambig', ...
        'c1',NaN,'errA',NaN,'errB',NaN,'dBIC',NaN,'f1',NaN,'f2',NaN,'tSep',NaN,'overlap',NaN, ...
        'atomThr',NaN,'included',false,'reason','');
end

function [P,res] = mpfit(xa, Q, C)
    P = []; res = [];
    try
        [P,res] = mp_adapt_chirplets(xa, Q, C.chirplet.M, C.chirplet.D, C.chirplet.i0, ...
            C.chirplet.radix, 'no', C.chirplet.mnits, C.chirplet.level, ...
            'RefineAlgorithm','expectmax', 'PType','Oneill');
    catch
    end
end

function [traj, nChain, nPar] = geometry(atoms, C)
    nSig = size(atoms,1); nChain=0; nPar=0; traj=NaN;
    if nSig < 2
        if nSig==1, traj = 1; end   % a lone atom: treat as fully-chained (geometry decided in classify)
        return;
    end
    medSig = median(atoms(:,4));
    for j = 1:nSig-1
        tcj=atoms(j,1); tck=atoms(j+1,1); fcj=atoms(j,2); fck=atoms(j+1,2);
        sj=atoms(j,4); sk=atoms(j+1,4); crj=atoms(j,3); crk=atoms(j+1,3);
        gap = tck - tcj;
        overlapSpan = min(tcj+sj, tck+sk) - max(tcj-sj, tck-sk);
        df = fck - fcj;                          % next minus current
        % parallel edge (dual-like): temporally COEXIST + distinct fc + both near-flat -> take priority
        if overlapSpan > 0 && abs(df) > C.chirplet.dualDf && abs(crj) <= C.chirplet.cTolHzs && abs(crk) <= C.chirplet.cTolHzs
            nPar = nPar + 1;
        % chain edge (single-like): temporally adjacent + fc NON-INCREASING (descend OR flat plateau,
        % so a steep-drop-then-plateau nonlinear sweep still reads as ONE connected trajectory).
        elseif gap <= C.chirplet.chainGapK*medSig && df <= C.chirplet.dfTol
            nChain = nChain + 1;
        end
    end
    if (nChain+nPar) > 0, traj = nChain/(nChain+nPar); else, traj = NaN; end
end

function cl = classify(atoms, nChain, nPar, C)
    nSig = size(atoms,1);
    if nSig == 0, cl='ambig'; return; end
    if nSig == 1
        if atoms(1,3) <= -C.chirplet.crMin, cl='single'; else, cl='ambig'; end  % lone flat atom = ambig
        return;
    end
    if nPar >= 1, cl='dual'; return; end
    if ~isnan(nChain) && (nChain/max(nChain+nPar,1)) >= C.chirplet.trajThiSingle && nPar==0
        cl='single'; return;
    end
    cl='ambig';
end

function [f1,f2,tSep,overlap] = descriptors(atoms, cl)
    f1=NaN; f2=NaN; tSep=NaN; overlap=NaN;
    if size(atoms,1) < 1, return; end
    if strcmp(cl,'dual') && size(atoms,1)>=2
        % the two highest-energy distinct-freq atoms
        [~,o] = sort(atoms(:,5),'descend'); a=atoms(o(1),:); b=atoms(o(2),:);
        f1 = max(a(2),b(2)); f2 = min(a(2),b(2)); tSep = abs(a(1)-b(1));
        overlap = min(a(1)+a(4), b(1)+b(4)) - max(a(1)-a(4), b(1)-b(4));
    else
        f1 = atoms(1,2); f2 = atoms(end,2); tSep = atoms(end,1)-atoms(1,1); overlap = NaN;
    end
end
