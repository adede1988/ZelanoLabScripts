function [y, fsOut, eventTable] = generateParametricWaveFromW( ...
    W, bandEdges, elementID, onsetSec, intensity, durationSec, fsOut, opts)
% generateParametricWaveFromW
%
% Generates a synthetic wave sound from NMF spectral templates.
%
% Inputs:
%   W           = frequency x component matrix from result.synth.W
%   bandEdges   = frequency band edges from result.synth.bandEdges
%   elementID   = which component each event uses
%   onsetSec    = onset time of each event, seconds
%   intensity   = relative POWER of each event
%   durationSec = duration of each event, seconds
%   fsOut       = output sample rate
%   opts        = optional parameter struct
%
% The function immediately plays the sound and does not save anything.
%
% Example:
%   [y, fs] = generateParametricWaveFromW( ...
%       result.synth.W, ...
%       result.synth.bandEdges, ...
%       [1 3 5 7], ...
%       [0 1.2 1.8 3.0], ...
%       [0.2 0.5 0.2 0.1], ...
%       [3.0 1.0 1.5 2.5], ...
%       12000);
%
% If elementID is empty, the function assumes one event per component:
%   elementID = 1:size(W,2)

if nargin < 7 || isempty(fsOut)
    fsOut = 12000;
end

if nargin < 8
    opts = struct;
end

% ----------------------------
% Defaults
% ----------------------------
opts = setDefault_local(opts, "masterDurationSec", []);
opts = setDefault_local(opts, "normalizeOutput", true);
opts = setDefault_local(opts, "targetPeak", 0.95);
opts = setDefault_local(opts, "playSound", true);

opts = setDefault_local(opts, "attackSec", 1);
opts = setDefault_local(opts, "releaseSec", 0.5);
opts = setDefault_local(opts, "envelopeType", "attackDecay");
opts = setDefault_local(opts, "decayCurve", 3.0);

opts = setDefault_local(opts, "stereo", false);
opts = setDefault_local(opts, "randomSeed", []);

opts = setDefault_local(opts, "lowFreqRollOffHz", 20);
opts = setDefault_local(opts, "highFreqRollOffHz", fsOut / 2);

opts = setDefault_local(opts, "backgroundNoiseFloor", 0.10);
% backgroundNoiseFloor can add a very small always-on shaped bed.
% Keep at 0 for pure event synthesis.

if ~isempty(opts.randomSeed)
    rng(opts.randomSeed);
end

% ----------------------------
% Validate inputs
% ----------------------------
nComponents = size(W, 2);

if isempty(elementID)
    elementID = 1:nComponents;
end

elementID = elementID(:);
onsetSec = onsetSec(:);
intensity = intensity(:);
durationSec = durationSec(:);

nEvents = numel(elementID);

if numel(onsetSec) ~= nEvents || ...
   numel(intensity) ~= nEvents || ...
   numel(durationSec) ~= nEvents

    error("elementID, onsetSec, intensity, and durationSec must have the same length.");
end

if any(elementID < 1) || any(elementID > nComponents)
    error("elementID contains values outside the range of W components.");
end

if any(durationSec <= 0)
    error("All durationSec values must be positive.");
end

if any(onsetSec < 0)
    error("All onsetSec values must be >= 0.");
end

% Intensities are interpreted as relative POWER and normalized to sum to 1
intensity(intensity < 0) = 0;

if sum(intensity) == 0
    intensity(:) = 1;
end

intensityPower = intensity ./ sum(intensity);

% Convert power weights to amplitude weights
intensityAmp = sqrt(intensityPower);

% ----------------------------
% Output duration
% ----------------------------
if isempty(opts.masterDurationSec)
    totalDurSec = max(onsetSec + durationSec) + 0.25;
else
    totalDurSec = opts.masterDurationSec;
end

nOut = ceil(totalDurSec * fsOut);
y = zeros(nOut, 1);

% ----------------------------
% Generate and sum events
% ----------------------------
for e = 1:nEvents

    k = elementID(e);

    nEvent = max(8, round(durationSec(e) * fsOut));

    event = synthesizeComponentNoise_local( ...
        W(:, k), ...
        bandEdges, ...
        nEvent, ...
        fsOut, ...
        opts);

    env = makeTemporalEnvelope_local(nEvent, fsOut, opts);
    event = event .* env;

    % Normalize event RMS before applying requested power weight
    eventRMS = sqrt(mean(event.^2)) + eps;
    event = event ./ eventRMS;

    event = event * intensityAmp(e);

    startSamp = round(onsetSec(e) * fsOut) + 1;
    stopSamp = min(nOut, startSamp + nEvent - 1);

    validN = stopSamp - startSamp + 1;

    if validN > 0
        y(startSamp:stopSamp) = y(startSamp:stopSamp) + event(1:validN);
    end
end

% ----------------------------
% Optional always-on background noise bed
% ----------------------------
if opts.backgroundNoiseFloor > 0

    meanSpectrum = mean(W, 2);
    bed = synthesizeComponentNoise_local( ...
        meanSpectrum, ...
        bandEdges, ...
        nOut, ...
        fsOut, ...
        opts);

    bed = bed ./ (sqrt(mean(bed.^2)) + eps);
    y = y + opts.backgroundNoiseFloor * bed;
end

% ----------------------------
% Final normalization
% ----------------------------
if opts.normalizeOutput
    peakVal = max(abs(y)) + eps;
    y = y ./ peakVal * opts.targetPeak;
end

% Optional stereo duplication
if opts.stereo
    y = [y y];
end

% ----------------------------
% Event table
% ----------------------------
eventTable = table;
eventTable.eventNumber = (1:nEvents)';
eventTable.elementID = elementID;
eventTable.onsetSec = onsetSec;
eventTable.durationSec = durationSec;
eventTable.offsetSec = onsetSec + durationSec;
eventTable.inputIntensity = intensity;
eventTable.normalizedPower = intensityPower;
eventTable.amplitudeMultiplier = intensityAmp;

% ----------------------------
% Play immediately
% ----------------------------
if opts.playSound
    sound(y, fsOut);
end

end


% ============================================================
% Create shaped noise from one component spectrum
% ============================================================
function x = synthesizeComponentNoise_local(componentSpectrum, bandEdges, nSamples, fs, opts)

componentSpectrum = componentSpectrum(:);
componentSpectrum(componentSpectrum < 0) = 0;

if max(componentSpectrum) > 0
    componentSpectrum = componentSpectrum ./ max(componentSpectrum);
else
    componentSpectrum(:) = 1;
end

bandCenters = sqrt(bandEdges(1:end-1) .* bandEdges(2:end));
bandCenters = bandCenters(:);

% Convert relative power-like template into amplitude template
ampTemplate = sqrt(componentSpectrum + eps);

% FFT length
nfft = 2^nextpow2(nSamples);

white = randn(nfft, 1);
X = fft(white);

% Frequency vector folded around Nyquist
freq = (0:nfft-1)' * fs / nfft;
foldedFreq = min(freq, fs - freq);

% Interpolate template onto FFT bins
ampResponse = interp1( ...
    bandCenters, ...
    ampTemplate, ...
    foldedFreq, ...
    "pchip", ...
    0);

% Frequency rolloffs
ampResponse(foldedFreq < opts.lowFreqRollOffHz) = 0;
ampResponse(foldedFreq > opts.highFreqRollOffHz) = 0;

% Avoid DC offset
ampResponse(1) = 0;

Y = X .* ampResponse;

x = real(ifft(Y));
x = x(1:nSamples);

x = x - mean(x);
x = x ./ (sqrt(mean(x.^2)) + eps);

end


% ============================================================
% Temporal envelope
% ============================================================
function env = makeTemporalEnvelope_local(nSamples, fs, opts)

durSec = nSamples / fs;
t = linspace(0, durSec, nSamples)';

attackSec = min(opts.attackSec, durSec * 0.45);
releaseSec = min(opts.releaseSec, durSec * 0.45);

attackN = max(1, round(attackSec * fs));
releaseN = max(1, round(releaseSec * fs));

switch string(opts.envelopeType)

    case "flat"
        env = ones(nSamples, 1);

    case "hann"
        env = hann(nSamples, "periodic");

    case "attackDecay"
        env = ones(nSamples, 1);

        % Smooth attack
        attack = linspace(0, 1, attackN)';
        attack = 0.5 - 0.5 * cos(pi * attack);
        env(1:attackN) = attack;

        % Exponential-ish decay over full duration
        decay = exp(-opts.decayCurve * t / max(durSec, eps));
        decay = decay ./ max(decay + eps);

        env = env .* decay;

        % Smooth release
        release = linspace(1, 0, releaseN)';
        release = 0.5 - 0.5 * cos(pi * release);
        env(end-releaseN+1:end) = env(end-releaseN+1:end) .* release;

    case "swellDecay"
        env = ones(nSamples, 1);

        swell = sin(pi * min(t / max(durSec * 0.35, eps), 1) / 2);
        decay = exp(-opts.decayCurve * max(t - durSec * 0.25, 0) / max(durSec, eps));

        env = swell .* decay;

        release = linspace(1, 0, releaseN)';
        release = 0.5 - 0.5 * cos(pi * release);
        env(end-releaseN+1:end) = env(end-releaseN+1:end) .* release;

    otherwise
        error("Unknown envelopeType: %s", opts.envelopeType);
end

env = env ./ max(env + eps);

end


% ============================================================
% Defaults helper
% ============================================================
function opts = setDefault_local(opts, fieldName, value)

if ~isfield(opts, fieldName) || isempty(opts.(fieldName))
    opts.(fieldName) = value;
end

end