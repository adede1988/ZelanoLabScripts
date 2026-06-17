function result = decomposeOceanWaveElements(wavFile, outDir)
% decomposeOceanWaveElements
%
% Decomposes a long ocean-wave recording into recurring acoustic elements
% using spectrogram-band features + nonnegative matrix factorization.
%
% The output is:
%   1. components_table.csv
%   2. element_events.csv
%   3. diagnostic figures
%
% Example:
%   result = decomposeOceanWaveElements( ...
%       "G:\My Drive\beach_waves_1hour.wav", ...
%       "G:\My Drive\wave_decomp_output");

if nargin < 2
    outDir = "wave_decomp_output";
end

if ~exist(outDir, "dir")
    mkdir(outDir);
end

% ----------------------------
% Parameters
% ----------------------------
targetFs = 12000;       % analysis sample rate
nElements = 10;         % number of acoustic elements to extract

winSec = 0.20;          % spectrogram window
hopSec = 0.05;          % spectrogram hop

nBands = 48;            % compressed frequency bands
fMin = 20;              % lowest frequency to analyze
minEventDurSec = 0.15;  % discard shorter activations
mergeGapSec = 0.20;     % merge nearby activations

activationZThresh = 1.5; % threshold for component activation
nmfIterations = 250;

% ----------------------------
% Load audio
% ----------------------------
fprintf("Loading audio...\n");

[x, fsOrig] = audioread(wavFile);

% Convert to mono
if size(x, 2) > 1
    x = mean(x, 2);
end

x = x - mean(x);
x = x ./ max(abs(x) + eps);

durationSec = numel(x) / fsOrig;

fprintf("Original sample rate: %.1f Hz\n", fsOrig);
fprintf("Duration: %.2f minutes\n", durationSec / 60);

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

win = hann(winN, "periodic");

[S, F, T] = spectrogram(xAnalysis, win, noverlap, nfft, fs);

P = abs(S).^2; % power spectrogram

% ----------------------------
% Compress spectrogram into log-spaced bands
% ----------------------------
fprintf("Compressing into frequency bands...\n");

fMax = fs / 2;
bandEdges = logspace(log10(fMin), log10(fMax), nBands + 1);

V = zeros(nBands, numel(T));

for b = 1:nBands
    idx = F >= bandEdges(b) & F < bandEdges(b + 1);

    if any(idx)
        V(b, :) = mean(P(idx, :), 1);
    end
end

% Remove all-zero bands if any
badBands = all(V == 0, 2);
V(badBands, :) = [];
bandEdgesForTable = bandEdges;
nBandsActual = size(V, 1);

% Log compression
V = log10(V + eps);

% Baseline subtract each band so that NMF focuses on changing structure
for b = 1:size(V, 1)
    V(b, :) = V(b, :) - prctile(V(b, :), 10);
end

% Ensure nonnegative
V(V < 0) = 0;

% Normalize each band so low and high frequencies both contribute
bandScale = median(V, 2) + eps;
V = V ./ bandScale;

% Light temporal smoothing
smoothFrames = max(3, round(0.15 / hopSec));
V = smoothdata(V, 2, "gaussian", smoothFrames);

% ----------------------------
% NMF decomposition
% ----------------------------
fprintf("Running NMF with %d elements...\n", nElements);

[W, H] = simpleNMF(V, nElements, nmfIterations);

% Normalize components so activations are comparable
for k = 1:nElements
    scale = max(W(:, k)) + eps;
    W(:, k) = W(:, k) / scale;
    H(k, :) = H(k, :) * scale;
end

% ----------------------------
% Characterize components
% ----------------------------
fprintf("Characterizing elements...\n");

componentInfo = table;

bandCenters = sqrt(bandEdges(1:end-1) .* bandEdges(2:end));
bandCenters = bandCenters(:);

% In case bad bands were removed
if numel(bandCenters) ~= size(W, 1)
    bandCenters = linspace(fMin, fMax, size(W, 1))';
end

for k = 1:nElements

    specProfile = W(:, k);
    act = H(k, :);

    spectralCentroid = sum(specProfile .* bandCenters) / sum(specProfile + eps);

    lowEnergy = sum(specProfile(bandCenters < 250));
    midEnergy = sum(specProfile(bandCenters >= 250 & bandCenters < 2000));
    highEnergy = sum(specProfile(bandCenters >= 2000));

    activationSparsity = mean(act > prctile(act, 90));

    suggestedLabel = suggestWaveElementLabel( ...
        spectralCentroid, lowEnergy, midEnergy, highEnergy, activationSparsity);

    componentInfo.elementID(k, 1) = k;
    componentInfo.suggestedLabel(k, 1) = string(suggestedLabel);
    componentInfo.spectralCentroidHz(k, 1) = spectralCentroid;
    componentInfo.lowBandWeight(k, 1) = lowEnergy;
    componentInfo.midBandWeight(k, 1) = midEnergy;
    componentInfo.highBandWeight(k, 1) = highEnergy;
    componentInfo.activationSparsity(k, 1) = activationSparsity;
end

writetable(componentInfo, fullfile(outDir, "components_table.csv"));

% ----------------------------
% Detect activation events for each element
% ----------------------------
fprintf("Detecting element activations...\n");

allEvents = table;

for k = 1:nElements

    act = H(k, :);

    % Robust z-score
    medAct = median(act);
    madAct = 1.4826 * median(abs(act - medAct)) + eps;
    actZ = (act - medAct) / madAct;

    % Smooth activation
    actZs = smoothdata(actZ, "gaussian", max(3, round(0.20 / hopSec)));

    active = actZs > activationZThresh;

    eventsK = activeMaskToEvents(active, T, actZs, k, ...
        componentInfo.suggestedLabel(k), minEventDurSec, mergeGapSec);

    allEvents = [allEvents; eventsK];
end

% Sort by onset
if ~isempty(allEvents)
    allEvents = sortrows(allEvents, "onsetSec");
end

writetable(allEvents, fullfile(outDir, "element_events.csv"));

% ----------------------------
% Diagnostic plots
% ----------------------------
fprintf("Making diagnostic plots...\n");

% Plot element spectral templates
fig1 = figure("Color", "w", "Position", [100 100 1200 700]);
tiledlayout(5, 2, "TileSpacing", "compact");

for k = 1:nElements
    nexttile;
    semilogx(bandCenters, W(:, k), "LineWidth", 2);
    grid on;
    xlabel("Frequency, Hz");
    ylabel("Weight");
    title(sprintf("Element %d: %s", k, componentInfo.suggestedLabel(k)));
end

saveas(fig1, fullfile(outDir, "element_spectral_profiles.png"));

% Plot activations over time
fig2 = figure("Color", "w", "Position", [100 100 1400 800]);
offset = 0;

hold on;

for k = 1:nElements
    act = H(k, :);
    act = act ./ max(act + eps);
    plot(T / 60, act + offset, "LineWidth", 1.2);
    text(T(end) / 60, offset + 0.4, ...
        sprintf("Element %d: %s", k, componentInfo.suggestedLabel(k)), ...
        "FontSize", 9);
    offset = offset + 1.2;
end

xlabel("Time, minutes");
ylabel("Normalized activation + offset");
title("Activation timeline for extracted wave elements");
grid on;

saveas(fig2, fullfile(outDir, "element_activation_timeline.png"));

% Optional: plot sequence of detected events
if ~isempty(allEvents)
    fig3 = figure("Color", "w", "Position", [100 100 1400 600]);
    hold on;

    for i = 1:height(allEvents)
        y = allEvents.elementID(i);
        plot([allEvents.onsetSec(i), allEvents.offsetSec(i)] / 60, ...
             [y, y], ...
             "LineWidth", 4);
    end

    yticks(1:nElements);
    yticklabels("E" + string(1:nElements));
    xlabel("Time, minutes");
    ylabel("Element");
    title("Detected sequence of acoustic wave elements");
    grid on;

    saveas(fig3, fullfile(outDir, "detected_element_sequence.png"));
end

% ----------------------------
% Return result struct
% ----------------------------
result.wavFile = wavFile;
result.outDir = outDir;
result.fsOriginal = fsOrig;
result.fsAnalysis = fs;
result.durationSec = durationSec;
result.componentInfo = componentInfo;
result.events = allEvents;
result.W = W;
result.H = H;
result.T = T;
result.bandCenters = bandCenters;

fprintf("Done.\n");
fprintf("Wrote output to: %s\n", outDir);

end


% ============================================================
% Simple NMF using multiplicative updates
% ============================================================
function [W, H] = simpleNMF(V, K, nIter)

rng(1);

[F, N] = size(V);

W = rand(F, K);
H = rand(K, N);

epsVal = 1e-9;

for iter = 1:nIter

    % Update H
    H = H .* ((W' * V) ./ (W' * W * H + epsVal));

    % Update W
    W = W .* ((V * H') ./ (W * (H * H') + epsVal));

    % Normalize occasionally
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
% Convert active mask to event table
% ============================================================
function events = activeMaskToEvents(active, T, actZ, elementID, label, minDurSec, mergeGapSec)

events = table;

active = active(:);
actZ = actZ(:);

d = diff([false; active; false]);
starts = find(d == 1);
stops = find(d == -1) - 1;

if isempty(starts)
    return
end

onsets = T(starts);
offsets = T(stops);

% Merge events separated by short gaps
mergedOnsets = [];
mergedOffsets = [];

curOn = onsets(1);
curOff = offsets(1);

for i = 2:numel(onsets)
    gap = onsets(i) - curOff;

    if gap <= mergeGapSec
        curOff = offsets(i);
    else
        mergedOnsets = [mergedOnsets; curOn];
        mergedOffsets = [mergedOffsets; curOff];

        curOn = onsets(i);
        curOff = offsets(i);
    end
end

mergedOnsets = [mergedOnsets; curOn];
mergedOffsets = [mergedOffsets; curOff];

% Build table
row = 0;

for i = 1:numel(mergedOnsets)

    dur = mergedOffsets(i) - mergedOnsets(i);

    if dur < minDurSec
        continue
    end

    idx = T >= mergedOnsets(i) & T <= mergedOffsets(i);

    if ~any(idx)
        continue
    end

    [peakZ, relIdx] = max(actZ(idx));
    idxList = find(idx);
    peakTime = T(idxList(relIdx));

    row = row + 1;

    events.elementID(row, 1) = elementID;
    events.suggestedLabel(row, 1) = string(label);
    events.onsetSec(row, 1) = mergedOnsets(i);
    events.peakSec(row, 1) = peakTime;
    events.offsetSec(row, 1) = mergedOffsets(i);
    events.durationSec(row, 1) = dur;
    events.peakActivationZ(row, 1) = peakZ;
end

end


% ============================================================
% Rough semantic label suggestion
% ============================================================
function label = suggestWaveElementLabel(centroidHz, lowEnergy, midEnergy, highEnergy, sparsity)

total = lowEnergy + midEnergy + highEnergy + eps;

lowFrac = lowEnergy / total;
midFrac = midEnergy / total;
highFrac = highEnergy / total;

if lowFrac > 0.60 && sparsity > 0.20
    label = "low background surf bed";
elseif lowFrac > 0.55 && sparsity <= 0.20
    label = "deep incoming rumble";
elseif midFrac > 0.45 && highFrac > 0.25 && centroidHz < 1800
    label = "broadband crash";
elseif highFrac > 0.50 && sparsity < 0.15
    label = "continuous foamy hiss";
elseif highFrac > 0.50 && sparsity >= 0.15
    label = "short splash or ripple";
elseif centroidHz < 300
    label = "low surge or wash";
elseif centroidHz < 1000
    label = "lapping water movement";
elseif centroidHz < 2500
    label = "whitewater decay";
else
    label = "high-frequency sparkle";
end

end