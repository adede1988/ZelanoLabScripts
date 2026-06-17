function result = decomposeOceanWaveElements_v2(wavFile, nElements)
% decomposeOceanWaveElements_v2
%
% Decomposes a long ocean recording into NMF spectral components and returns
% the W matrix in a form ready for synthetic sound generation.
%
% Main output for synthesis:
%   result.synth.W
%   result.synth.bandEdges
%   result.synth.bandCenters
%   result.synth.fsRecommended
%
% Example:
%   result = decomposeOceanWaveElements_v2("beach_waves.wav", 10);

if nargin < 2
    nElements = 10;
end

% ----------------------------
% Analysis parameters
% ----------------------------
targetFs = 12000;

winSec = 0.20;
hopSec = 0.05;

nBands = 64;
fMin = 20;
fMax = targetFs / 2;

nmfIterations = 300;

% ----------------------------
% Load audio
% ----------------------------
fprintf("Loading audio...\n");

[x, fsOrig] = audioread(wavFile);

if size(x, 2) > 1
    x = mean(x, 2);
end

x = x - mean(x);
x = x ./ max(abs(x) + eps);

durationSec = numel(x) / fsOrig;

fprintf("Original Fs: %.1f Hz\n", fsOrig);
fprintf("Duration: %.2f min\n", durationSec / 60);

% Downsample for analysis
if fsOrig ~= targetFs
    xAnalysis = resample(x, targetFs, fsOrig);
else
    xAnalysis = x;
end

fs = targetFs;

% ----------------------------
% Spectrogram
% ----------------------------
fprintf("Computing spectrogram...\n");

winN = round(winSec * fs);
hopN = round(hopSec * fs);
noverlap = winN - hopN;
nfft = 2^nextpow2(winN);

[S, F, T] = spectrogram( ...
    xAnalysis, ...
    hann(winN, "periodic"), ...
    noverlap, ...
    nfft, ...
    fs);

P = abs(S).^2;

% ----------------------------
% Compress into log-spaced frequency bands
% ----------------------------
fprintf("Compressing spectrogram into bands...\n");

bandEdges = logspace(log10(fMin), log10(fMax), nBands + 1);
bandCenters = sqrt(bandEdges(1:end-1) .* bandEdges(2:end));

V = zeros(nBands, numel(T));

for b = 1:nBands
    idx = F >= bandEdges(b) & F < bandEdges(b + 1);

    if any(idx)
        V(b, :) = mean(P(idx, :), 1);
    end
end

% Log compression
Vlog = log10(V + eps);

% Baseline subtract each frequency band
for b = 1:nBands
    Vlog(b, :) = Vlog(b, :) - prctile(Vlog(b, :), 10);
end

% Make nonnegative for NMF
Vlog(Vlog < 0) = 0;

% Normalize each frequency band so high-frequency components are not ignored
bandScale = median(Vlog, 2) + eps;
Vnmf = Vlog ./ bandScale;

% Smooth over time
smoothFrames = max(3, round(0.15 / hopSec));
Vnmf = smoothdata(Vnmf, 2, "gaussian", smoothFrames);

% ----------------------------
% NMF
% ----------------------------
fprintf("Running NMF with %d elements...\n", nElements);

[W, H] = simpleNMF_local(Vnmf, nElements, nmfIterations);

% Normalize W/H pairs
for k = 1:nElements
    s = max(W(:, k)) + eps;
    W(:, k) = W(:, k) / s;
    H(k, :) = H(k, :) * s;
end

% Make a synthesis-friendly W matrix
% Rows = frequency bands
% Columns = components
Wsynth = W;

for k = 1:nElements
    Wsynth(:, k) = Wsynth(:, k) ./ max(Wsynth(:, k) + eps);
end

% ----------------------------
% Component summary
% ----------------------------
componentInfo = table;

for k = 1:nElements

    spec = Wsynth(:, k);
    centroidHz = sum(spec(:) .* bandCenters(:)) / sum(spec(:) + eps);

    lowWeight = sum(spec(bandCenters < 250));
    midWeight = sum(spec(bandCenters >= 250 & bandCenters < 2000));
    highWeight = sum(spec(bandCenters >= 2000));

    componentInfo.elementID(k, 1) = k;
    componentInfo.spectralCentroidHz(k, 1) = centroidHz;
    componentInfo.lowWeight(k, 1) = lowWeight;
    componentInfo.midWeight(k, 1) = midWeight;
    componentInfo.highWeight(k, 1) = highWeight;
    componentInfo.suggestedLabel(k, 1) = string( ...
        suggestWaveElementLabel_local(centroidHz, lowWeight, midWeight, highWeight));
end

% ----------------------------
% Diagnostic plots
% ----------------------------
figure("Color", "w", "Position", [100 100 1200 750]);
tiledlayout(ceil(nElements / 2), 2, "TileSpacing", "compact");

for k = 1:nElements
    nexttile;
    semilogx(bandCenters, Wsynth(:, k), "LineWidth", 2);
    grid on;
    xlabel("Frequency, Hz");
    ylabel("Relative weight");
    title(sprintf("Element %d: %s", k, componentInfo.suggestedLabel(k)));
end

figure("Color", "w", "Position", [100 100 1400 750]);
hold on;

offset = 0;

for k = 1:nElements
    act = H(k, :);
    act = act ./ max(act + eps);

    plot(T / 60, act + offset, "LineWidth", 1);
    text(T(end) / 60, offset + 0.4, ...
        sprintf("E%d: %s", k, componentInfo.suggestedLabel(k)), ...
        "FontSize", 9);

    offset = offset + 1.2;
end

xlabel("Time, min");
ylabel("Normalized activation + offset");
title("NMF component activations");
grid on;

% ----------------------------
% Output
% ----------------------------
result = struct;

result.wavFile = wavFile;
result.fsOriginal = fsOrig;
result.fsAnalysis = fs;
result.durationSec = durationSec;

result.W = W;
result.H = H;
result.T = T;

result.Vnmf = Vnmf;

result.bandEdges = bandEdges;
result.bandCenters = bandCenters;

result.componentInfo = componentInfo;

% This is the main structure for synthesis
result.synth = struct;
result.synth.W = Wsynth;
result.synth.bandEdges = bandEdges;
result.synth.bandCenters = bandCenters;
result.synth.fsRecommended = targetFs;
result.synth.componentInfo = componentInfo;

fprintf("Done.\n");
fprintf("Use result.synth.W and result.synth.bandEdges for sound generation.\n");

end


% ============================================================
% Local NMF helper
% ============================================================
function [W, H] = simpleNMF_local(V, K, nIter)

rng(1);

[nFreq, nTime] = size(V);

W = rand(nFreq, K);
H = rand(K, nTime);

epsVal = 1e-9;

for iter = 1:nIter

    H = H .* ((W' * V) ./ (W' * W * H + epsVal));
    W = W .* ((V * H') ./ (W * (H * H') + epsVal));

    if mod(iter, 25) == 0
        for k = 1:K
            s = sum(W(:, k)) + epsVal;
            W(:, k) = W(:, k) / s;
            H(k, :) = H(k, :) * s;
        end

        fprintf("  NMF iteration %d / %d\n", iter, nIter);
    end
end

end


% ============================================================
% Label helper
% ============================================================
function label = suggestWaveElementLabel_local(centroidHz, lowWeight, midWeight, highWeight)

total = lowWeight + midWeight + highWeight + eps;

lowFrac = lowWeight / total;
midFrac = midWeight / total;
highFrac = highWeight / total;

if lowFrac > 0.60
    label = "low surf bed or deep rumble";
elseif lowFrac > 0.40 && midFrac > 0.35
    label = "surge or incoming wash";
elseif midFrac > 0.50 && highFrac < 0.30
    label = "broad crash body";
elseif midFrac > 0.35 && highFrac > 0.35
    label = "crash plus splash";
elseif highFrac > 0.60 && centroidHz > 2500
    label = "foam hiss or sparkle";
elseif highFrac > 0.45
    label = "ripples or lapping";
else
    label = "mixed wave texture";
end

end