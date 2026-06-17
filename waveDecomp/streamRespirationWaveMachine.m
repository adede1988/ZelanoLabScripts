function [eventTable, streamState] = streamRespirationWaveMachine(respTrace, respFs, W, bandEdges, opts)
% streamRespirationWaveMachine
%
% Real-time / closed-loop style respiration-controlled wave synthesizer.
%
% This does NOT pre-render the full audio file.
% It walks through the respiration trace at 60 Hz and streams audio block-by-block.
%
% Inputs:
%   respTrace : respiration timeseries
%   respFs    : respiration sampling rate
%   W         : NMF component spectra, frequency x component
%   bandEdges : frequency band edges corresponding to rows of W
%   opts      : optional settings struct
%
% Example:
%   [eventTable, streamState] = streamRespirationWaveMachine( ...
%       respTrace, respFs, result.W, result.bandEdges);

if nargin < 5
    opts = struct;
end

% ------------------------------------------------------------
% Defaults
% ------------------------------------------------------------
opts = setDefault(opts, "controlFs", 60);
opts = setDefault(opts, "fsOut", 12000);

opts = setDefault(opts, "historySec", 30);
opts = setDefault(opts, "grainDurSec", 5.0);
opts = setDefault(opts, "atackSec", .5);
opts = setDefault(opts, "releaseSec", .2);
opts = setDefault(opts, "decayCurve", 1.0);




opts = setDefault(opts, "baseTargetPeak", 0.20);
opts = setDefault(opts, "minTargetPeak", 0.0005);
opts = setDefault(opts, "sameBinDecay", 0.96);
opts = setDefault(opts, "newBinPenalty", 0.04);

opts = setDefault(opts, "nGrainVariants", 16);
opts = setDefault(opts, "plotWindowSec", 20);
opts = setDefault(opts, "plotUpdateHz", 15);

opts = setDefault(opts, "centerRespiration", true);
opts = setDefault(opts, "softClip", true);
opts = setDefault(opts, "limiterDrive", 1.5);

opts = setDefault(opts, "randomSeed", []);

if ~isempty(opts.randomSeed)
    rng(opts.randomSeed);
end

% ------------------------------------------------------------
% Check for streaming audio support
% ------------------------------------------------------------
if exist("audioDeviceWriter", "class") ~= 8 && exist("audioDeviceWriter", "file") ~= 2
    error(["audioDeviceWriter was not found. This streaming version requires " + ...
           "Audio Toolbox / DSP-style real-time audio output. Base MATLAB sound() " + ...
           "cannot do low-latency closed-loop streaming cleanly."]);
end

% ------------------------------------------------------------
% Fixed bin-to-component mapping
% ------------------------------------------------------------
% binCodes     = [-3 -2 -1  1  2  3  4  5  6];
  binCodes     = [-4 -3 -2 -1  1  2  3  4  5];
  componentMap = [ 9  2  4  5  10 6  7  1  3]; %low freq exhale
% componentMap = [ 4 10  2  5  8  3  7  9  6]; %low freq inhale

% ------------------------------------------------------------
% Prepare respiration as a simulated live input stream
% ------------------------------------------------------------
respTrace = respTrace(:);
respTrace = fillmissing(respTrace, "linear", "EndValues", "nearest");

if opts.centerRespiration
    respTrace = respTrace - median(respTrace, "omitnan");
end

tResp = (0:numel(respTrace)-1)' / respFs;
totalDurSec = tResp(end);

controlFs = opts.controlFs;
tControl = (0:1/controlFs:totalDurSec)';
nTicks = numel(tControl);

% This is only for simulating the live stream from a stored trace.
% In real hardware use, replace this with your DAQ sample read.
respControl = interp1(tResp, respTrace, tControl, "linear", "extrap");

% ------------------------------------------------------------
% Audio streaming setup
% ------------------------------------------------------------
fsOut = opts.fsOut;
audioBlockN = round(fsOut / controlFs);

if abs(audioBlockN - fsOut / controlFs) > 1e-9
    error("fsOut must be divisible by controlFs. For example, 12000 / 60 = 200 samples per tick.");
end

deviceWriter = audioDeviceWriter( ...
    "SampleRate", fsOut, ...
    "BufferSize", audioBlockN);

grainN = round(opts.grainDurSec * fsOut);

% Rolling mix buffer. New grains are added at the front; each loop writes
% the first audioBlockN samples, then shifts the buffer left.
mixBufferN = grainN + 4 * audioBlockN;
mixBuffer = zeros(mixBufferN, 1);

% ------------------------------------------------------------
% Precompute component grains
% ------------------------------------------------------------
fprintf("Precomputing component grains...\n");

nComponents = size(W, 2);
grainBank = cell(nComponents, opts.nGrainVariants);

for k = 1:nComponents
    for v = 1:opts.nGrainVariants
        bb = find(componentMap == k); 
        %custom longer grains for pause time: 
        if binCodes(bb) == -1 | binCodes(bb) == 1
            tmp = opts; 
            tmp.grainDurSec = tmp.grainDurSec +3; 
        else
            tmp = opts; 
        end

        grainBank{k, v} = synthesizeComponentGrain( ...
            W(:, k), ...
            bandEdges, ...
            grainN, ...
            fsOut, ...
            tmp);
    end
end

% ------------------------------------------------------------
% Live figure
% ------------------------------------------------------------
fig = figure("Color", "w", "Position", [100 100 1200 600]);
ax = axes(fig);
hold(ax, "on");

hResp = plot(ax, NaN, NaN, "k", "LineWidth", 1.5);
hNow = xline(ax, 0, "r", "LineWidth", 2);

% Bin-limit horizontal lines
nBinEdges = 10;  % -3/-2/-1/0/1/2/3/4/5/6 boundaries
hBinLimits = gobjects(nBinEdges, 1);

for b = 1:nBinEdges
    hBinLimits(b) = yline(ax, NaN, ":", ...
        "LineWidth", 1, ...
        "Alpha", 0.45);
end

yline(ax, 0, "--", "Color", [0.5 0.5 0.5]);

xlabel(ax, "Time, sec");
ylabel(ax, "Respiration");
title(ax, "Streaming respiration-controlled wave synthesis");
grid(ax, "off");

hText = text(ax, 0.02, 0.95, "", ...
    "Units", "normalized", ...
    "VerticalAlignment", "top", ...
    "FontSize", 12, ...
    "FontWeight", "bold", ...
    "BackgroundColor", "w");

% ------------------------------------------------------------
% Event logging
% ------------------------------------------------------------
eventLog = struct;
eventLog.timeSec = nan(nTicks, 1);
eventLog.respValue = nan(nTicks, 1);
eventLog.localMax = nan(nTicks, 1);
eventLog.localMin = nan(nTicks, 1);
eventLog.binCode = nan(nTicks, 1);
eventLog.componentID = nan(nTicks, 1);
eventLog.repeatCount = nan(nTicks, 1);
eventLog.previousRepeatCount = nan(nTicks, 1);
eventLog.targetPeak = nan(nTicks, 1);
eventLog.variantID = nan(nTicks, 1);

% ------------------------------------------------------------
% Streaming state
% ------------------------------------------------------------
lastBin = NaN;
repeatCount = 0;
previousRepeatCount = 0;
%avoid double playing components due to quick jitter: 
countBack = zeros(30,length(binCodes)); 




historyN = round(opts.historySec * controlFs);
plotUpdateEvery = max(1, round(controlFs / opts.plotUpdateHz));

fprintf("Starting streaming playback...\n");

% Prime audio device
setup(deviceWriter, zeros(audioBlockN, 1));

try
    for i = 1:nTicks

        tNow = tControl(i);
        rNow = respControl(i);
 
        % ----------------------------------------------------
        % Use previous 30 sec of respiration to get local range
        % ----------------------------------------------------
        histStart = max(1, i - historyN + 1);
        histVals = respControl(histStart:i);

        localMax = prctile(histVals, 97); %max(histVals);
        localMin = prctile(histVals, 10); %min(histVals);

        if localMax <= 0
            localMax = max(abs(histVals)) + eps;
        end

        if localMin >= 0
            localMin = -max(abs(histVals)) - eps;
        end

        
        % ----------------------------------------------------
        % Bin current respiration value
        % ----------------------------------------------------
        binNow = respirationValueToBin(rNow, localMax, localMin);
        
        mapIdx = find(binCodes == binNow, 1);
        componentNow = componentMap(mapIdx);

        % ----------------------------------------------------
        % Repeat-dependent intensity scaling
        % ----------------------------------------------------
        if isnan(lastBin)

            repeatCount = 1;
            previousRepeatCount = 0;
            targetPeak = opts.baseTargetPeak;

        elseif binNow == lastBin

            repeatCount = repeatCount + 1;

            % Same bin/component: weaker with each repeat.
            targetPeak = targetPeak * opts.sameBinDecay^(repeatCount - 1);
            if binNow == 1 || binNow == -1
                targetPeak = targetPeak * .1; 
            end
            

        else

            if binNow - lastBin > 1 %jump by more than one bin
                binNow = binNow - 1; 
                if binNow == 0
                    binNow = 1; 
                end
            end

            if binNow - lastBin < -1 %jump by more than one bin
                binNow = binNow + 1; 
                if binNow == 0
                    binNow = -1; 
                end
            end

            previousRepeatCount = repeatCount;
            repeatCount = 1;

            % New bin/component: weakened according to how many repeats
            % happened in the previous component.

            targetPeak = opts.baseTargetPeak / ...
                (1 + opts.newBinPenalty * previousRepeatCount);

            if lastBin>0 && binNow < lastBin && binNow > 0 && binNow ~= 1%heading down the inhale
                targetPeak = 0; 
            end

            if lastBin<0 && binNow > lastBin && binNow < 0 && binNow ~= -1%heading up the exhale
                targetPeak = 0; 
            end

            

        end

        countBack(1:29, :) = countBack(2:30, :); 
        countBack(30, :) = binCodes == binNow; 

        targetPeak = max(opts.minTargetPeak, min(1, targetPeak));
        if (binNow == -1 || binNow == 1 ) & ... % if you're going to play pause sound
                sum(countBack(:, binCodes==-1 | binCodes==1) == 1, 'all') > 20 % and you already have a lot
                targetPeak = 0;  %then don't! 
        end

        if (binNow ~= -1 & binNow ~= 1 ) & ...% if any other sound
                sum(countBack(:, binCodes==binNow)) > 5 %just don't play it too too much
                targetPeak = 0; 
        end

        % ----------------------------------------------------
        % Create one new 1-second component grain
        % ----------------------------------------------------
        variantID = randi(opts.nGrainVariants);
        grain = grainBank{componentNow, variantID};

        grain = grain ./ (max(abs(grain)) + eps);
        grain = grain * targetPeak;

        % Add new grain into rolling mix buffer at current time.
        mixBuffer(1:grainN) = mixBuffer(1:grainN) + grain;

        % ----------------------------------------------------
        % Output the next short audio block immediately
        % ----------------------------------------------------
        audioBlock = mixBuffer(1:audioBlockN);

        if opts.softClip
            audioBlock = softClipAudio(audioBlock, opts.limiterDrive);
        end

        deviceWriter(audioBlock);

        % Advance rolling buffer
        mixBuffer(1:end-audioBlockN) = mixBuffer(audioBlockN+1:end);
        mixBuffer(end-audioBlockN+1:end) = 0;

        % ----------------------------------------------------
        % Log event
        % ----------------------------------------------------
        eventLog.timeSec(i) = tNow;
        eventLog.respValue(i) = rNow;
        eventLog.localMax(i) = localMax;
        eventLog.localMin(i) = localMin;
        eventLog.binCode(i) = binNow;
        eventLog.componentID(i) = componentNow;
        eventLog.repeatCount(i) = repeatCount;
        eventLog.previousRepeatCount(i) = previousRepeatCount;
        eventLog.targetPeak(i) = targetPeak;
        eventLog.variantID(i) = variantID;

        lastBin = binNow;

        % ----------------------------------------------------
        % Live visualization
        % ----------------------------------------------------
        if mod(i, plotUpdateEvery) == 0 || i == 1

            winStart = max(0, tNow - opts.plotWindowSec);
            plotIdx = tControl >= winStart & tControl <= tNow;
            
            set(hResp, ...
                "XData", tControl(plotIdx), ...
                "YData", respControl(plotIdx));
            
            hNow.Value = tNow;
            
            xlim(ax, [winStart, max(winStart + opts.plotWindowSec, tNow + 0.1)]);
            
            % ------------------------------------------------------------
            % Plot current adaptive bin limits
            % ------------------------------------------------------------
            negEdges = linspace(localMin, 0, 4);     % localMin, -2, -1, 0 boundaries
            posEdges = linspace(0, localMax, 5);     % 0, 1, 2, 3, 4, 5, 6 boundaries
            
            % Combine without duplicating zero
            binEdges = [negEdges(1:end-1), 0, posEdges(2:end)];
            
            for b = 1:numel(binEdges)
                hBinLimits(b).Value = binEdges(b);
            end
            
            % ------------------------------------------------------------
            % Y limits should include both trace and adaptive bin edges
            % ------------------------------------------------------------
            yMin = min([respControl(plotIdx); binEdges(:)]);
            yMax = max([respControl(plotIdx); binEdges(:)]);
            
            if yMin == yMax
                yMin = yMin - 1;
                yMax = yMax + 1;
            end
            
            pad = 0.1 * (yMax - yMin);
            ylim(ax, [yMin - pad, yMax + pad]);
            
            hText.String = sprintf( ...
                "t = %.2f s\nresp = %.3f\nbin = %d\ncomponent = %d\nrepeat = %d\ntargetPeak = %.3f\nlocalMin = %.3f\nlocalMax = %.3f", ...
                tNow, ...
                rNow, ...
                binNow, ...
                componentNow, ...
                repeatCount, ...
                targetPeak, ...
                localMin, ...
                localMax);
            
            drawnow limitrate nocallbacks;
        end

        if ~isvalid(fig)
            break
        end
    end

catch ME
    release(deviceWriter);
    rethrow(ME);
end

release(deviceWriter);

% ------------------------------------------------------------
% Convert log to table
% ------------------------------------------------------------
eventTable = struct2table(eventLog);
eventTable = eventTable(~isnan(eventTable.timeSec), :);

streamState = struct;
streamState.opts = opts;
streamState.binCodes = binCodes;
streamState.componentMap = componentMap;
streamState.fsOut = fsOut;
streamState.audioBlockN = audioBlockN;
streamState.grainN = grainN;

fprintf("Streaming finished.\n");

end

% ========================================================================
% Helper: assign current respiration value to bin
% ========================================================================
function binNow = respirationValueToBin(rNow, localMax, localMin)

if rNow >= 0

    % Positive space from 0 to max split into 6 bins:
    % 1, 2, 3, 4, 5
    if localMax <= 0
        binNow = 1;
        return
    end

    edges = linspace(0, localMax, 5);
    binNow = discretize(rNow, edges);

    if isnan(binNow)
        if rNow >= localMax
            binNow = 5;
        else
            binNow = 1;
        end
    end

    binNow = max(1, min(5, binNow));

else

    % Negative space from 0 to min split into 3 bins:
    % -1, -2, -3, -4
    maxNegMag = abs(localMin);

    if maxNegMag <= 0
        binNow = -1;
        return
    end

    mag = abs(rNow);
    edges = linspace(0, maxNegMag, 5);

    negLevel = discretize(mag, edges);

    if isnan(negLevel)
        if mag >= maxNegMag
            negLevel = 4;
        else
            negLevel = 1;
        end
    end

    negLevel = max(1, min(4, negLevel));
    binNow = -negLevel;
end

end

% ========================================================================
% Helper: synthesize one component grain from W spectrum
% ========================================================================
function grain = synthesizeComponentGrain(componentSpectrum, bandEdges, nSamples, fs, opts)

componentSpectrum = componentSpectrum(:);
componentSpectrum(componentSpectrum < 0) = 0;

if max(componentSpectrum) > 0
    componentSpectrum = componentSpectrum ./ max(componentSpectrum);
else
    componentSpectrum(:) = 1;
end

bandCenters = sqrt(bandEdges(1:end-1) .* bandEdges(2:end));
bandCenters = bandCenters(:);

% Treat W as power-like; convert to amplitude weighting.
ampTemplate = sqrt(componentSpectrum + eps);

nfft = 2^nextpow2(nSamples);

white = randn(nfft, 1);
X = fft(white);

freq = (0:nfft-1)' * fs / nfft;
foldedFreq = min(freq, fs - freq);

ampResponse = interp1( ...
    bandCenters, ...
    ampTemplate, ...
    foldedFreq, ...
    "pchip", ...
    0);

% Basic cleanup
ampResponse(foldedFreq < 20) = 0;
ampResponse(foldedFreq > fs/2) = 0;
ampResponse(1) = 0;

Y = X .* ampResponse;

grain = real(ifft(Y));
grain = grain(1:nSamples);

grain = grain - mean(grain);
grain = grain ./ (sqrt(mean(grain.^2)) + eps);

% Envelope: fast attack, exponential decay, soft release.
t = linspace(0, 1, nSamples)';

attackSec = opts.attackSec; %.5
releaseSec = opts.releaseSec; % .2
decayCurve = opts.decayCurve; %2.5

attackN = max(1, round(attackSec * fs));
releaseN = max(1, round(releaseSec * fs));

env = exp(-decayCurve * t);
env = env ./ max(env + eps);

attack = linspace(0, 1, attackN)';
attack = 0.5 - 0.5 * cos(pi * attack);
env(1:attackN) = env(1:attackN) .* attack;

release = linspace(1, 0, releaseN)';
release = 0.5 - 0.5 * cos(pi * release);
env(end-releaseN+1:end) = env(end-releaseN+1:end) .* release;

grain = grain .* env;
grain = grain ./ (max(abs(grain)) + eps);

end

% ========================================================================
% Helper: soft clipping for safe streaming
% ========================================================================
function y = softClipAudio(x, drive)

y = tanh(drive * x) ./ tanh(drive);

end

% ========================================================================
% Helper: defaults
% ========================================================================
function opts = setDefault(opts, fieldName, value)

if ~isfield(opts, fieldName) || isempty(opts.(fieldName))
    opts.(fieldName) = value;
end

end