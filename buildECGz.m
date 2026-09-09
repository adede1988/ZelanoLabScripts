function [ECGz, beatSep] = buildECGz(outDat)
%BUILDECGZ  Band-pass + z-score the ECG channels and return the beat-separation.
%
%   [ECGz, beatSep] = buildECGz(outDat)
%
%   Shared by processECG (beat detection -> HRV) and paramCheckECG (interactive
%   verification of the beat-detection spec) so both operate on identical ECGz
%   and beatSep. Breathing task only -- it is the only task with ECG channels.
%
%   ECGz    : channels whose label contains 'ECG', band-passed 5-40 Hz and
%             z-scored per channel (rows = channels, cols = time).
%   beatSep : minimum separation (samples) enforced between detected beats.

    idx = cellfun(@(x) contains(x, 'ECG'), outDat.labels);
    ECG = outDat.data(idx, :);

    d = designfilt('bandpassiir', 'FilterOrder', 4, ...
        'HalfPowerFrequency1', 5, 'HalfPowerFrequency2', 40, ...
        'SampleRate', outDat.fs);
    ECG = filtfilt(d, ECG')';

    ECGz = (ECG - mean(ECG, 2)) ./ std(ECG, [], 2);

    beatSep = outDat.fs / 20;
end
