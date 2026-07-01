function [M, relTimeMs] = v2_peakcenter(mat, tMs, gpt, halfMs, fs)
% V2_PEAKCENTER  Re-center each trial's time series to its own gammaPeakTime.
%   [M, relTimeMs] = v2_peakcenter(mat, tMs, gpt, halfMs, fs)
%   mat : [nTrial x nT] time series on the tMs grid (finalOnset-locked)
%   gpt : [nTrial x 1] gammaPeakTime (ms rel finalOnset); NaN -> NaN row
%   Returns M [nTrial x nPk] on a relative axis +-halfMs (0 = gammaPeakTime), NaN off-window.

    nTr = size(mat,1);
    half = round(halfMs/1000*fs); nPk = 2*half+1;
    relTimeMs = (-half:half)/fs*1000;
    M = nan(nTr, nPk);
    for i = 1:nTr
        if ~isfinite(gpt(i)), continue; end
        [~,k] = min(abs(tMs - gpt(i)));
        lo = k-half; hi = k+half;
        srcLo = max(1,lo); srcHi = min(numel(tMs),hi);
        dstLo = srcLo-lo+1; dstHi = dstLo+(srcHi-srcLo);
        M(i,dstLo:dstHi) = mat(i, srcLo:srcHi);
    end
end
