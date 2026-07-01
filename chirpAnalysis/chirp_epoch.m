function E = chirp_epoch(sig, fs, onsets, C)
% CHIRP_EPOCH  Build padded sniff-locked epochs around finalOnset for one channel.
%   E = chirp_epoch(sig, fs, onsets, C)
%
%   sig    : 1 x nSamp continuous signal (one channel, e.g. bestMac or a macBP)
%   onsets : trial onset samples (behDat.finalOnset), @ fs
%   C      : chirp_config (uses C.epochWin and C.padSec)
%
%   Returns epochs WITH padding (C.padSec each side) so the caller can FIR-filter / Hilbert /
%   TFR without edge artifacts, then trim to the core window via E.coreIdx. Trials whose
%   padded window runs off the recording are marked ~E.valid (kept as NaN rows so trial
%   indexing stays aligned with behDat rows).
%
%   Fields:
%     E.dataPad [nTrial x nPad]  padded epochs (NaN rows where invalid)
%     E.tPadMs  [1 x nPad]       time (ms) of padded samples rel. onset
%     E.coreIdx [1 x nCore] logical into padded cols spanning C.epochWin
%     E.tMs     [1 x nCore]      time (ms) of core samples rel. onset
%     E.valid   [nTrial x 1] logical
%     E.onsets  [nTrial x 1] (echoed)  E.fs

    if nargin < 4 || isempty(C), C = chirp_config(); end
    sig = sig(:)';  T = numel(sig);  N = numel(onsets);

    s0p = round((C.epochWin(1) - C.padSec) * fs);   % padded start offset (samples, <=0)
    s1p = round((C.epochWin(2) + C.padSec) * fs);   % padded end   offset
    s0c = round(C.epochWin(1) * fs);                % core start offset
    s1c = round(C.epochWin(2) * fs);                % core end   offset
    nPad = s1p - s0p + 1;

    dataPad = nan(N, nPad);
    valid   = false(N, 1);
    for i = 1:N
        e = onsets(i);
        if ~isfinite(e) || e <= 0, continue; end
        a = round(e) + s0p;  b = round(e) + s1p;
        if a < 1 || b > T, continue; end
        seg = sig(a:b);
        if any(~isfinite(seg)), continue; end
        dataPad(i, :) = seg;
        valid(i) = true;
    end

    padSampOffsets = s0p:s1p;                 % offsets of padded cols rel. onset
    coreIdx = padSampOffsets >= s0c & padSampOffsets <= s1c;

    E = struct();
    E.dataPad = dataPad;
    E.tPadMs  = padSampOffsets / fs * 1000;
    E.coreIdx = coreIdx;
    E.tMs     = padSampOffsets(coreIdx) / fs * 1000;
    E.valid   = valid;
    E.onsets  = onsets(:);
    E.fs      = fs;
    E.padSec  = C.padSec;
    E.epochWin= C.epochWin;
end
