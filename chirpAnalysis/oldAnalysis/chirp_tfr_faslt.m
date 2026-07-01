function tfr = chirp_tfr_faslt(E, C)
% CHIRP_TFR_FASLT  FASLT power TFR for all trials of one channel (cached kernels).
%   tfr = chirp_tfr_faslt(E, C)
%   E : epoch struct from chirp_epoch (E.dataPad [nTrial x nPad], E.coreIdx, E.tMs, E.valid, E.fs)
%   C : chirp_config
%
%   Builds the FASLT bank ONCE (chirp_faslt_bank), applies it per valid trial on the PADDED
%   epoch, then trims to the [-1,+3]s core. Output power tensor + trial-mean (power-averaged,
%   preserves induced gamma). The ONLY licensed use of FASLT power is here (spectrogram display
%   + ridge). Phase/envelope for the other tests come from hilbert() of the FIR-bandpassed signal.
%
%   tfr fields: .power [Nf x nCore x nTrial] (NaN slabs for invalid trials)
%               .meanPower [Nf x nCore]  .freqs [1 x Nf] .tMs [1 x nCore] .fs
%               .bandMask (freqs in C.ridge.band) .valid .bank (cached) .padOK

    fs = E.fs;
    bank = chirp_faslt_bank(fs, C);

    % F1 pad assertion: epoch pad must cover the max wavelet half-support in the RIDGE band
    fb = bank.F >= C.ridge.band(1) & bank.F <= C.ridge.band(2);
    halfSupp = 0;
    for i = find(fb)
        for k = 1:bank.order_int(i)
            if ~isempty(bank.wavelets{i,k}), halfSupp = max(halfSupp, fix(numel(bank.wavelets{i,k})/2)); end
        end
    end
    padSamp = round(C.padSec*fs);
    tfr.padOK = padSamp >= halfSupp;
    if ~tfr.padOK
        warning('chirp_tfr_faslt:pad', ...
            'padSec=%.2fs (%d samp) < max ridge-band wavelet half-support %d samp (%.2fs). Raise C.padSec.', ...
            C.padSec, padSamp, halfSupp, halfSupp/fs);
    end

    [nTrial, nPad] = size(E.dataPad);
    Nf = numel(bank.F);
    nCore = sum(E.coreIdx);
    P = nan(Nf, nCore, nTrial);
    for i = 1:nTrial
        if ~E.valid(i), continue; end
        wt = chirp_faslt_apply(E.dataPad(i,:), bank);   % [Nf x nPad] power
        P(:,:,i) = wt(:, E.coreIdx);
    end

    tfr.power    = P;
    tfr.meanPower= mean(P(:,:,E.valid), 3, 'omitnan');
    tfr.freqs    = bank.F;
    tfr.tMs      = E.tMs;
    tfr.fs       = fs;
    tfr.bandMask = bank.F >= C.ridge.band(1) & bank.F <= C.ridge.band(2);
    tfr.valid    = E.valid;
    tfr.bank     = bank;
end
