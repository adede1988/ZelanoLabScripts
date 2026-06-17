clear

%% Paths + session list (carried over from cueTask loader)

codePre   = 'C:\Users\Adam\Documents\GitHub\';

% Centralized behavioral results location for "newSet" participants
behDatPath_newSet = 'R:\Neurology\Zelano_Lab\Lab_Common\OBE_task_backup\tasks\pea_threshold\results';






 

%% Toolboxes / paths (kept from original PEA script; adjust as needed on your machine)
addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
addpath(genpath([codePre 'myFrequentUse']))
addpath([codePre 'myFrequentUse/export_fig_repo'])
addpath(genpath('C:\Users\Adam\Documents\eeglab2026.0.0'))

addpath([codePre 'fieldtrip-20230118'])
addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])

set(0, 'defaultfigurewindowstyle', 'docked')
% ft_defaults  % FieldTrip not installed on this machine and not used by this script

%% Main loop
cfg        = applyParams('threshTask','makeOutDat');
sessionIDs = cfg.sessionIDs;
datPre     = cfg.datPre;
datPrei    = cfg.datPrei;
newSet     = cfg.newIDs;
rspIDX     = cfg.rspIDX;
rspFlip    = cfg.rspFlip;

for sessi = 1:numel(sessionIDs)
try
if ~isempty(getenv('RUNSUBSET')) && ~ismember(sessi, str2num(getenv('RUNSUBSET'))), continue; end

    sessID = sessionIDs{sessi};
    disp(['Working on: ' sessID])

  preProcDir = fullfile(datPre{datPrei(sessi)}, sessID, 'preProc');


 if true || exist(fullfile(preProcDir, [sessID '_PEA_threshold_preproc.mat']), 'file')
    %% -----------------------
    %  Load RAW data (.mat)
    %% -----------------------
    rawRoot = fullfile(datPre{datPrei(sessi)}, sessID, 'raw');
    datFolders = dir(rawRoot);
    datFolders = datFolders([datFolders.isdir]);

    idx = find(cellfun(@(x) contains(x, 'raw_PEAintensityPleasantness'), {datFolders.name}));

    if ~isscalar(idx)
        error('PEA loader: expected exactly 1 raw_PEAintensityPleasantness folder for %s, found %d.', ...
            sessID, numel(idx));
    end

    rawMatPath = fullfile(datFolders(idx).folder, datFolders(idx).name, [datFolders(idx).name '.mat']);
    if ~exist(rawMatPath, 'file')
        error('PEA loader: raw .mat not found for %s at %s', sessID, rawMatPath);
    end

    dat = load(rawMatPath);
    if ~isfield(dat, 'curDat')
        error('PEA loader: raw file for %s does not contain curDat.', sessID);
    end
    dat = dat.curDat;

    %% -----------------------
    %  Load BEHAVIORAL data
    %  - newSet: centralized results folder (behDatPath_newSet\<sessID>\*.mat)
    %  - old set: search within participant's local "behavior" folder
    %% -----------------------
    isNew = ismember(sessID, newSet);

    if isNew
        behDir = dir(fullfile(behDatPath_newSet, sessID));
        behDir = behDir(~[behDir.isdir]);
        matIdx = find(contains({behDir.name}, '.mat'));

        if ~isscalar(matIdx)
            error('PEA loader: expected exactly 1 behavioral .mat for newSet %s in %s, found %d.', ...
                sessID, fullfile(behDatPath_newSet, sessID), numel(matIdx));
        end

        behFile = fullfile(behDir(matIdx).folder, behDir(matIdx).name);

    else
        % Find behavioral folder inside subject directory (folder name containing 'ehavior')
        subFolders = dir(fullfile(datPre{datPrei(sessi)}, sessID));
        subFolders = subFolders([subFolders.isdir]);
        % remove '.' and '..'
        subFolders = subFolders(~ismember({subFolders.name}, {'.','..'}));

        behFolderIdx = find(contains(lower({subFolders.name}), 'ehavior'));
        if ~isscalar(behFolderIdx)
            error('PEA loader: behavioral folder not uniquely identified for %s under %s (found %d matches).', ...
                sessID, fullfile(datPre{datPrei(sessi)}, sessID), numel(behFolderIdx));
        end

        behRoot = fullfile(subFolders(behFolderIdx).folder, subFolders(behFolderIdx).name);

        % Recursively search for mat files that look like PEA threshold outputs
        behMats = dir(fullfile(behRoot, '**', '*.mat'));
        if isempty(behMats)
            error('PEA loader: no behavioral .mat files found for %s under %s', sessID, behRoot);
        end

        % Filter to likely PEA-related files (tighten this filter if needed)
        keep = false(size(behMats));
        for k = 1:numel(behMats)
            nm = lower(behMats(k).name);
            pth = lower(behMats(k).folder);
            keep(k) = contains(nm, 'pea') || contains(nm, 'pleasant') || contains(nm, 'intensity') || contains(nm, 'threshold') || ...
                      contains(pth, 'pea') || contains(pth, 'pleasant') || contains(pth, 'intensity') || contains(pth, 'threshold');
        end
        behMats = behMats(keep);

        if ~isscalar(behMats)
            % Provide a helpful error listing candidates
            msg = sprintf('PEA loader: expected exactly 1 PEA behavioral .mat for %s under %s, found %d.\nCandidates:\n', ...
                sessID, behRoot, numel(behMats));
            for k = 1:numel(behMats)
                msg = [msg sprintf('  - %s\n', fullfile(behMats(k).folder, behMats(k).name))]; %#ok<AGROW>
            end
            error('%s', msg);
        end

        behFile = fullfile(behMats.folder, behMats.name);
    end

    % Load + format behavioral table
    behLoaded = load(behFile);
    if ~isfield(behLoaded, 'outMat')
        error('PEA loader: behavioral file for %s does not contain outMat: %s', sessID, behFile);
    end

    behDat = behLoaded.outMat;
    behDat = cat(1, behDat{:});

    % Universal fill for missing pleasantness/intensity: 735
    f = @(~) {735};
    emptyIDX = cellfun(@isempty, behDat(:,4));
    behDat(emptyIDX,4) = cellfun(@(x) f(x), behDat(emptyIDX, 4));
    emptyIDX = cellfun(@isempty, behDat(:,5));
    behDat(emptyIDX,5) = cellfun(@(x) f(x), behDat(emptyIDX, 5));

    behDat = cell2mat(behDat);
    behDat = array2table(behDat, 'VariableNames', ...
        {'trialNum', 'Odor', 'RT', 'pleasantness', 'intensity', 'ITI', ...
         'Trial_StartTime', 'Trial_EndTime'});

    %% -----------------------
    %  Build outDat
    %% -----------------------
    outDat = struct;
    outDat.labels  = dat.outLabs;
    outDat.CSClist = dat.ncslabels;
    outDat.fs      = dat.rawData.fsample;
    outDat.sessID  = sessID;
    outDat.behDat  = behDat;

    % Keep these session-level params in case downstream code wants them
    outDat.rspIDX  = rspIDX(sessi);
    outDat.rspFlip = rspFlip(sessi);

    %% -----------------------
    %  TTL extraction (single event type; must be exactly 45)
    %% -----------------------
    idxEvt = cellfun(@(x) contains(x, 'event'), outDat.labels);
    if sum(idxEvt) ~= 1
        error('PEA loader: expected exactly 1 photodiode/event channel for %s, found %d.', ...
            sessID, sum(idxEvt));
    end

    photoDiode = dat.rawData.trial{1}(idxEvt, :);

    % Handle NaNs (photodiode + entire raw matrix)
    nanidx = isnan(photoDiode);
    nNan   = sum(nanidx);

    if nNan > 4000
        error('PEA loader: too many missing photodiode samples for %s (%d NaNs).', sessID, nNan);
    end

    rawData = dat.rawData.trial{1};

    if nNan > 0
        photoDiode = fillmissing(photoDiode, 'linear');
        photoDiode = fillmissing(photoDiode, 'nearest');

        rawData = fillmissing(rawData, 'linear', 2);
        rawData = fillmissing(rawData, 'nearest', 2);
    end

    outDat.data = rawData;

    % Normalize photodiode
    photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);

    % Subject-specific diode threshold (kept from original script)
    switch sessID
        case '250623_DUPI_NMH_KS_2'
            diodeThresh = -2;
        otherwise
            diodeThresh = -1.5;
    end

    downs = find(photoDiode(1:end-1) > diodeThresh & photoDiode(2:end) < diodeThresh);
    ups   = find(photoDiode(1:end-1) < diodeThresh & photoDiode(2:end) > diodeThresh);
    difVals = ups - downs;

    % Filter implausible pulse widths
    downs(difVals > 1500) = [];
    ups(difVals > 1500)   = [];
    difVals(difVals > 1500) = [];

    downs(difVals < 300) = [];
    ups(difVals < 300)   = [];
    difVals(difVals < 300) = [];

    % Trim off pre-start pulses (train detection)
    starti = 1;
    for di = 3:length(downs)
        if downs(di) - downs(di-2) < 3500
            starti = di;
        end
    end
    downs(1:starti)   = [];
    difVals(1:starti) = [];

    if length(downs) ~= 45
        error('PEA loader: wrong number of TTLs for %s (expected 45, found %d).', sessID, length(downs));
    end

    % Quick diagnostic plot
    figure
    plot(photoDiode)
    xline(downs, 'color', 'magenta', 'linewidth', 2)
    title(sprintf('Photodiode TTLs: %s (n=%d)', sessID, length(downs)))

    TTLs = round(downs(:) ./ 4);  % sample->ms index (assuming 4 kHz -> 1 kHz)

    %% -----------------------
    %  Metadata + save
    %% -----------------------
    outDat.task = "PEAintensityPleasantness_threshold";
    outDat.OGdataDir = fullfile(datPre{datPrei(sessi)}, sessID);

    tmp = dir(fullfile(datPre{datPrei(sessi)}, sessID, '*.m'));
    tmp = tmp(contains({tmp.name}, 'LoadData'));
    if numel(tmp) == 1
        outDat.loadFile = tmp.name;
    else
        tmp2 = dir(fullfile(datPre{datPrei(sessi)}, sessID, '*AD.m'));
        if numel(tmp2) == 1
            outDat.loadFile = tmp2.name;
        else
            error('PEA loader: load file not identified uniquely for %s', sessID);
        end
    end

    outDat.preProcScript = 'threshPreProc_makeOutDat.m';
    if datPrei(sessi) == 1
        outDat.type = 'Dupi';
    elseif datPrei(sessi) == 2
        outDat.type = 'OBE';
    else
        error('PEA loader: unexpected datPrei value for %s', sessID);
    end

    outDat.TTL = table;
    outDat.TTL.sniff = TTLs;

    preProcDir = fullfile(datPre{datPrei(sessi)}, sessID, 'preProc');
    if ~exist(preProcDir, 'dir')
        mkdir(preProcDir);
    end

    save(fullfile(preProcDir, [sessID '_PEA_threshold_preproc.mat']), 'outDat', '-v7.3');

    % Optional: clear big per-session variables
    clear datFolders idx rawMatPath dat behDir behLoaded behDat behFile ...
          photoDiode rawData downs ups difVals TTLs tmp tmp2 idxEvt
 end
catch MErun
    warning('RUNFAIL %s: %s', sessID, MErun.message);
end
clear outDat dat dat1 dat2 comboDat photoDiode rawData behDat behDat1 behDat2 TTLs TTLs2 savedTTL data
end
