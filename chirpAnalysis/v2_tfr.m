function T = v2_tfr(sig, fs, onsets, good, C, endLimit)
% V2_TFR  Epoch bestMac -> broadband analytic + FASLT TFR -> baseline z-score (spec 4).
%   T = v2_tfr(sig, fs, onsets, good, C [,endLimit])
%   endLimit (optional, Nx1 absolute sample): O15 long-window guard -- blank each epoch beyond
%   the next finalOnset (post-truncation TFR/padData -> NaN; burstTruncated flag set).
%   sig    : 1 x nSamp bestMac channel (continuous)
%   onsets : behDat.finalOnset samples
%   good   : N x 1 logical (valid & not-noisy) trials to process; others -> NaN
%
%   Steps: padded epochs (C.padSec each side) -> per-trial demean+detrend -> FASLT POWER
%   (cached bank, reused across trials) -> trim pad -> per-frequency bootstrap baseline
%   z-score (myChanZscore, baseline C.baselineMs). Also returns the padded demeaned data so
%   the phase stage can narrowband-filter with padding then trim.
%
%   T fields: padData[N x nPad] coreIdx tMs tPadMs freqs[1xNf] zTFR[Nf x nCore x N]
%             good valid nCore N fs

    E = chirp_epoch(sig, fs, onsets, C);          % padded epochs (C.epochWin, C.padSec)
    N = numel(onsets); nPad = size(E.dataPad,2);
    coreIdx = E.coreIdx; tMs = E.tMs; nCore = numel(tMs);
    good = good(:) & E.valid;                       % only on-recording, non-noisy trials

    % per-trial demean + detrend (padded)
    padData = nan(N, nPad);
    for i = find(good)'
        x = E.dataPad(i,:);
        padData(i,:) = detrend(x - mean(x,'omitnan'));
    end

    % FASLT bank (built once, reused) -> per-trial POWER, trimmed to core
    bank = chirp_faslt_bank(fs, C);
    Nf = numel(bank.F);
    tfrPow = nan(Nf, nCore, N);
    for i = find(good)'
        wt = chirp_faslt_apply(padData(i,:), bank);   % [Nf x nPad] POWER
        tfrPow(:,:,i) = wt(:, coreIdx);
    end

    % per-frequency bootstrap baseline z-score (myChanZscore); baseline C.baselineMs
    bIdx = find(tMs >= C.baselineMs(1) & tMs <= C.baselineMs(2));
    basePeriod = [bIdx(1) bIdx(end)];
    gi = find(good)';
    zTFR = nan(Nf, nCore, N);
    if ~isempty(gi)
        for f = 1:Nf
            dataF = reshape(tfrPow(f, :, gi), nCore, numel(gi));  % [time x goodTrials]
            zF = myChanZscore(dataF, basePeriod);                 % [time x goodTrials]
            zTFR(f, :, gi) = reshape(zF, 1, nCore, numel(gi));
        end
    end

    % O15 long-window guard: blank each epoch beyond the next finalOnset
    burstTruncated = false(N,1);
    if nargin>=6 && ~isempty(endLimit)
        for i = find(good)'
            if ~isfinite(endLimit(i)), continue; end
            cutMs = (endLimit(i)-onsets(i))/fs*1000;
            zc = tMs >= cutMs;         if any(zc), zTFR(:, zc, i) = NaN; end
            pc = E.tPadMs >= cutMs;    if any(pc), padData(i, pc) = NaN; end
            burstTruncated(i) = cutMs <= C.epochWin(2)*1000;
        end
    end

    T = struct('padData',padData,'coreIdx',coreIdx,'tMs',tMs,'tPadMs',E.tPadMs, ...
        'freqs',bank.F(:)','zTFR',zTFR,'good',good,'valid',E.valid,'nCore',nCore, ...
        'N',N,'fs',fs,'burstTruncated',burstTruncated);
end
