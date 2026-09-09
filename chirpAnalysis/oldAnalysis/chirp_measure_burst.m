function out = chirp_measure_burst(S, C)
% CHIRP_MEASURE_BURST  Real-data burst-length + in-band SNR for one session (D4 + Phase-0 calib).
%   out = chirp_measure_burst(S, C)   % S from chirp_load_session
%
%   Gamma power = abs(hilbert(25-58 FIR))^2 of bestMac, finalOnset-locked, core-trimmed. Burst
%   start/stop via chirp_burstlen (75th-pctile-of-cross-trial-peak threshold; search (0,3000] ms).
%   In-band SNR (dB) = 10*log10(median peak gamma power / median pre-onset baseline power) -> feeds
%   the Phase-0 synthesizer so its regime matches reality.
%
%   out: .sessID .B (chirp_burstlen struct) .snrDb .gp .tMs .nTrial

    fs = S.fs;
    E = chirp_epoch(S.bestSig, fs, S.onsets, C);
    bb = chirp_bbfilt(E.dataPad, fs, C.bbBand, C);
    gpPad = abs(hilbert(bb')').^2;          % [nTrial x nPad] gamma power
    gp = gpPad(:, E.coreIdx);
    B = chirp_burstlen(gp, E.tMs, C);

    tMs = E.tMs;
    blMask = tMs >= -700 & tMs <= -100;
    post   = tMs > 0 & tMs <= C.burst.searchMs(2);
    snrTr = nan(size(gp,1),1);
    for i = 1:size(gp,1)
        if ~E.valid(i), continue; end
        bl = median(gp(i,blMask),'omitnan'); pk = max(gp(i,post));
        if bl>0 && isfinite(pk), snrTr(i) = 10*log10(pk/bl); end
    end

    out = struct('sessID',S.sessID,'B',B,'snrDb',median(snrTr,'omitnan'), ...
        'gp',gp,'tMs',tMs,'nTrial',sum(E.valid),'snrTrial',snrTr);
end
