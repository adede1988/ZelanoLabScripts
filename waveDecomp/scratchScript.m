result = decomposeOceanWaveElements_v2( ...
    "C:\Users\dtf8829\Downloads\OCEAN WAVES 1H - MALLORCA - QUIET MUSIC\OCEAN WAVES 1H - MALLORCA - QUIET MUSIC.wav", ...
     10);

% result = decomposeOceanWaveElements_v2( ...
%     "E:\oceanSounds\freetousesounds - Royalty Free Ocean Sounds & Seaside Ambience - 42 AMBSea-LR_Japan-Ambience, Seaside Tajima Island Gravel Stones.mp3", ...
%      10);
outDir = "G:\My Drive\GitHub\experiment_EEGsync_v2";

% Read plain numeric CSVs, no headers
W = readmatrix(fullfile(outDir, "wave_W_new.csv"));
bandEdges = readmatrix(fullfile(outDir, "wave_bandEdges.csv"));
load('R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\250815_EEG_NWU_PP\preProc\250815_EEG_NWU_PP_breathingPreProc.mat')

rspDat = chanDat.data(36,:) * -1; 
opts = struct; 
opts.attackSec = .1; 
opts.releaseSec= .5; 
opts.grainDurSec= 1; 
opts.decayCurve = 1; 
[eventTable, streamState] = streamRespirationWaveMachine(rspDat(1:60000), 500, W, bandEdges, opts);


% ------------------------------------------------------------
% Save NMF wave-synthesis matrices for PsychoPy
% ------------------------------------------------------------



if ~exist(outDir, "dir")
    mkdir(outDir);
end

% Use result.W directly
W_to_save = result.W;
bandEdges_to_save = result.bandEdges(:);

% Optional: if using v3 and you want synthesis-normalized W instead:
% W_to_save = result.synth.W;
% bandEdges_to_save = result.synth.bandEdges(:);

% Save plain numeric CSVs, no headers
writematrix(W_to_save, fullfile(outDir, "wave_W.csv"));
writematrix(bandEdges_to_save, fullfile(outDir, "wave_bandEdges.csv"));

% Also save component labels/info if available
if isfield(result, "componentInfo")
    writetable(result.componentInfo, fullfile(outDir, "wave_componentInfo.csv"));
end

fprintf("Saved PsychoPy export files to:\n%s\n", outDir);


















result = decomposeOceanWaveLibrary( ...
    "C:\Users\dtf8829\Downloads\OCEAN WAVES 1H - MALLORCA - QUIET MUSIC\OCEAN WAVES 1H - MALLORCA - QUIET MUSIC.wav", ...
     10);


for ii = 1:5
    tic
opts = struct; 
opts.targetPeak = .2; 
[y, fsOut, eventTable] = generateParametricWaveFromW( ...
    result.W, result.bandEdges, [1,], ... % element
                                [0], ... % onset
                                [1], ... % intensity
                                [60], ... % duration
                                12000, ...
                                opts);    
[y, fsOut, eventTable] = generateParametricWaveFromW( ...
    result.W, result.bandEdges, [5,  8,   3,   7,   9,   6,   4, 10,  2, 10], ...% element
                                [1,  1.1, 1.2, 1.3, 1.4, 1.6, 3, 3.5, 4, 6], ... % onset
                                [1,  1,   1.5, 2,   2,   4,   2, 1,   1, 1], ... % intensity
                                [4,  4,  4,   4,   4,   2.5,  1, 1,   1, 1], ... % duration
                                12000);


% Reproducible randomness
rng('shuffle')   % or rng(1) if you want the same result every time

% ------------------------------------------------------------
% Fixed opening texture
% ------------------------------------------------------------
baseElements   = [5,  8,   3,   7,   9,   6];
baseOnsets     = [1,  1.1, 1.2, 1.3, 1.4, 1.6];
baseIntensity  = [1,  1,   1.5, 2,   2,   4];
baseDuration   = [4,  4,   4,   4,   4,   2.5];

% ------------------------------------------------------------
% Random alternating sequence from 3 to 8 sec
% ------------------------------------------------------------
possibleElements = [4 10 2];

tStart = 3;
tEnd   = 8;

randElements  = [];
randOnsets    = [];
randIntensity = [];
randDuration  = [];

t = tStart;
lastElement = NaN;

while t <= tEnd

    % Pick randomly from 4, 10, and 2, avoiding immediate repeats
    availableElements = possibleElements(possibleElements ~= lastElement);
    thisElement = availableElements(randi(numel(availableElements)));

    randElements(end+1) = thisElement;
    randOnsets(end+1) = t;

    % Random intensity between 0.5 and 1.5
    randIntensity(end+1) = 0.5 + rand * (1.5 - 0.5);

    % Random duration between 1 and 2.5 sec
    randDuration(end+1) = 1.0 + rand * (2.5 - 1.0);

    % Next onset step between 0.1 and 0.75 sec
    t = t + 0.1 + rand * (0.75 - 0.1);

    lastElement = thisElement;
end

% ------------------------------------------------------------
% Combine fixed and random parts
% ------------------------------------------------------------
elementID = [baseElements, randElements];
onsetSec = [baseOnsets, randOnsets];
intensity = [baseIntensity, randIntensity];
durationSec = [baseDuration, randDuration];

% Optional: inspect generated event sequence
eventDesign = table( ...
    elementID(:), ...
    onsetSec(:), ...
    intensity(:), ...
    durationSec(:), ...
    'VariableNames', {'elementID', 'onsetSec', 'intensity', 'durationSec'});

disp(eventDesign)

% ------------------------------------------------------------
% Generate and play sound
% ------------------------------------------------------------
[y, fsOut, eventTable] = generateParametricWaveFromW( ...
    result.W, ...
    result.bandEdges, ...
    elementID, ...
    onsetSec, ...
    intensity, ...
    durationSec, ...
    12000);


pause(12)
toc
end




trans = analyzeWaveElementTransitions(result, ...
    "G:\My Drive\wave_decomposition");