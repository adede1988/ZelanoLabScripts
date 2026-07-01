function [x, gt] = chirp_synth_trial(family, P, C)
% CHIRP_SYNTH_TRIAL  Synthesize one Phase-0 trial (spec 5) on the padded epoch grid.
%   [x, gt] = chirp_synth_trial(family, P, C)
%   family : 'singleLin' | 'singleNL' | 'dualIndep' | 'dualLocked'
%   P      : overlapMs, Df, ampRatio([hi lo] or scalar), fLo, snrDb, burstMs, fastDropMs, seed
%            [optional] jitterMs (per-trial timing jitter; default 60). For single the WHOLE burst
%            jitters together; for dualIndep each tone jitters independently (the decoupling the
%            temporal test 6.9 detects); dualLocked jitters the pair together.
%
%   Returns x [1 x nPad] on the SAME padded grid chirp_epoch produces (onset at t=0) + ground
%   truth gt. Gamma embedded in pink+white noise at the requested in-band (25-58 Hz) SNR.

    fs = C.fs;
    s0p = round((C.epochWin(1)-C.padSec)*fs); s1p = round((C.epochWin(2)+C.padSec)*fs);
    off = s0p:s1p; tMs = off/fs*1000; t = off/fs;       % t=0 at onset
    N = numel(t);
    rng(P.seed, 'twister');
    if ~isfield(P,'jitterMs') || isempty(P.jitterMs), jit = 0.060; else, jit = P.jitterMs/1000; end

    fLo = P.fLo; fHi = fLo + P.Df;
    if isscalar(P.ampRatio) && ~isnan(P.ampRatio), aHi=P.ampRatio; aLo=1;
    elseif isscalar(P.ampRatio), aHi=1; aLo=1;
    else, aHi=P.ampRatio(1); aLo=P.ampRatio(2); end
    burstS = P.burstMs/1000;
    trialJit = jit*randn;                               % shared (whole-burst) timing jitter

    sig = zeros(1, N);
    gt = struct('family',family,'isDual',false,'fHi',fHi,'fLo',fLo,'ampHi',aHi,'ampLo',aLo, ...
        'coexistWinMs',[NaN NaN],'transitionWinMs',[NaN NaN],'nGammaCyclesOverlap',NaN,'ifTrue',[]);

    switch family
        case 'singleLin'
            tau = t - trialJit; m = tau>=0 & tau<=burstS;
            IF = nan(1,N); IF(m) = fHi + (fLo-fHi).*(tau(m)/burstS);
            env = onoff(t, trialJit, burstS+trialJit, fs);
            ph = 2*pi*cumIF(IF, fs, m);
            sig = aHi.*env.*cos(ph + 2*pi*rand);
            gt.transitionWinMs = [trialJit burstS+trialJit]*1000; gt.coexistWinMs = gt.transitionWinMs;
            gt.ifTrue = IF;

        case 'singleNL'
            fd = P.fastDropMs/1000; tau = t - trialJit; m = tau>=0 & tau<=burstS;
            IF = nan(1,N); md = tau>=0 & tau<fd; mp = tau>=fd & tau<=burstS;
            IF(md) = fLo + (fHi-fLo).*exp(-3*tau(md)/fd);
            IF(mp) = fLo + (fHi-fLo).*exp(-3);
            env = onoff(t, trialJit, burstS+trialJit, fs);
            ph = 2*pi*cumIF(IF, fs, m);
            sig = aHi.*env.*cos(ph + 2*pi*rand);
            gt.transitionWinMs = [trialJit fd+trialJit]*1000; gt.coexistWinMs = [trialJit burstS+trialJit]*1000;
            gt.ifTrue = IF;

        case {'dualIndep','dualLocked'}
            gt.isDual = true; ov = P.overlapMs/1000; cen = burstS/2;
            if strcmp(family,'dualIndep'), jH=jit*randn; jL=jit*randn; else, jH=trialJit; jL=trialJit; end
            hiOn = [0+jH,        cen+ov/2+jH];
            loOn = [cen-ov/2+jL, burstS+jL];
            envHi = onoff(t, hiOn(1), hiOn(2), fs);
            envLo = onoff(t, loOn(1), loOn(2), fs);
            phHi = 2*pi*fHi.*t + 2*pi*rand;
            if strcmp(family,'dualLocked'), phLo = 2*pi*fLo.*t + 0.7; else, phLo = 2*pi*fLo.*t + 2*pi*rand; end
            sig = aHi.*envHi.*cos(phHi) + aLo.*envLo.*cos(phLo);
            cw = [max(hiOn(1),loOn(1)) min(hiOn(2),loOn(2))];   % coexistence = both ON
            gt.coexistWinMs = cw*1000; gt.transitionWinMs = cw*1000;
            gt.nGammaCyclesOverlap = max(0,diff(cw))*mean([fHi fLo]);
            gt.ifTrue = nan(1,N);
        otherwise
            error('chirp_synth_trial:family','unknown family %s', family);
    end

    x = addInbandNoise(sig, fs, C.bbBand, P.snrDb);
    gt.tMs = tMs; gt.fs = fs;
end

% ---- helpers ----
function e = onoff(t, t0, t1, fs) %#ok<INUSD>
    tau = 0.080; k = 4/tau;
    e = 1./(1+exp(-k*(t-t0))) .* 1./(1+exp(k*(t-t1)));
    e(t < t0-0.2 | t > t1+0.2) = 0;
end
function ph = cumIF(IF, fs, m)
    inst = IF; inst(~m) = 0; inst(isnan(inst)) = 0;
    ph = cumsum(inst)/fs;
end
function y = addInbandNoise(sig, fs, band, snrDb)
    N = numel(sig); w = randn(1,N);
    F = fft(w); fr = (0:N-1)/N*fs; fr(fr>fs/2) = fs - fr(fr>fs/2); fr(fr<1) = 1;
    pink = real(ifft(F ./ sqrt(fr)));
    noise = 0.6*pink/std(pink) + 0.4*w/std(w);
    bs = bandvar(sig, fs, band); bn = bandvar(noise, fs, band);
    if bn <= 0, bn = eps; end
    noise = noise * sqrt((bs/(10^(snrDb/10))) / bn);
    y = sig + noise;
end
function v = bandvar(x, fs, band)
    n = numel(x); X = fft(x); fr = (0:n-1)/n*fs;
    m = (fr>=band(1) & fr<=band(2)) | (fr>=fs-band(2) & fr<=fs-band(1));
    v = sum(abs(X(m)).^2)/n^2;
end
