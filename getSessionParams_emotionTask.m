function [raw, P] = getSessionParams_emotionTask(S)
% LOAD_SESSION  Load raw data & behavior for a session and aggregate
%               subject-specific parameters into P.
%
% Usage:
%   S.id   = '250904_OBE_NWU_TI';
%   S.root = 'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\';
%   [raw, P] = load_session(S);
%
% Outputs:
%   raw: struct with fields
%       .sessID      - session ID string
%       .fs_raw      - original sampling rate (Hz)
%       .data        - raw data matrix [nChan x nSamples]
%       .labels      - channel labels (cellstr)
%       .ncslabels   - NCS labels (if present)
%       .beh         - behavior table for this session
%       .paths.*     - loaded file paths
%       .type        - 'Dupi' or 'OBE' (heuristic)
%   P: struct of parameters (subject-specific + defaults)
%       .task, .type, .fs_target
%       .debug
%       .rspIDX, .rspFlip
%       .respThresh, .cuedBackBuff, .adjWin
%       .hasEEG, .spikeThresh, .spikeWin, .spikeClean
%       .pd.* (photodiode pulse parsing thresholds)
%       .ttl.* (TTL expectations & overrides)

    arguments
        S struct
    end

    %% --- Build paths & load raw ---
    sessID = char(S.id);
    root   = char(S.root);

    matPath = fullfile(root, sessID, 'raw', 'raw_EmotionalMovieTask', 'raw_EmotionalMovieTask.mat');
    if ~exist(matPath, 'file')
        error('load_session:MissingMat', 'Raw MAT not found: %s', matPath);
    end
    dat = load(matPath);
    if isfield(dat, 'curDat')
        dat = dat.curDat;
    else
        error('load_session:BadMat', 'Expected curDat in MAT file: %s', matPath);
    end

   %no behavioral data needed for this task

   

    % Populate raw struct
    raw = struct();
    raw.sessID    = sessID;
    raw.fs_raw    = dat.rawData.fsample;
    raw.data      = dat.rawData.trial{1};

    raw.data(:, isnan(raw.data(end,:))) = []; 


    raw.labels    = dat.outLabs;
    if isfield(dat, 'ncslabels'), raw.ncslabels = dat.ncslabels; end
    raw.paths.mat = matPath;
    raw.paths.fig = fullfile(S.figPath, S.id);

    % Heuristic type (kept for bookkeeping)
    if contains(sessID, 'OBE', 'IgnoreCase', true) || contains(root, 'OBE', 'IgnoreCase', true)
        raw.type = 'OBE';
    else
        raw.type = 'Dupi';
    end

    %% --- Defaults (session-agnostic) ---
    P = struct();
    P.task      = 'EmotionalMovieTask';
    P.type      = raw.type;
    P.fs_target = 500;          % target downsample Fs
    P.debug     = false;        % turn plots on/off upstream

    % Resp channel selection (from your arrays: all 1's, all no-flip)
    P.rspIDX    = 1;
    P.rspFlip   = 1;


    % EEG / spike cleaning defaults
    P.hasEEG      = false;
    P.spikeThresh = 20;
    P.spikeWin    = 11;
    P.spikeClean  = false;

    % Photodiode pulse parsing 
    P.pd = struct();
    P.pd.zthresh        = -2;    % z-score threshold for "low"
    P.pd.minPulseSamp   = 350;   % remove pulses < 200 samples
    P.pd.maxPulseSamp   = 2000;  % remove pulses > 1200 samples
    P.pd.searchWin = 2000;
    P.pd.numNeg = 3;  
    P.pd.numPos = 2; 
    P.pd.numNeu = 1; 
    P.macroRemove = []; 
  

    %% --- Session-specific overrides (aggregated from your script) ---
    switch sessID
        % EEG present + spike cleaning parameters
        case '250904_OBE_NWU_TI'
            P.hasEEG      = true;
            P.spikeThresh = 50;
            P.spikeWin    = 11;
            P.spikeClean  = true;
            P.pd.zthresh  = -2;    % z-score threshold for "low"


           

        case '251009_OBE_NWU_CP_1'
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.spikeThresh = 15;
            P.spikeWin    = 7;
            P.macroRemove = [6]; 


        case '250225_OBE_NWU_AS_4'

        otherwise
            % keep defaults
    end

    % sanity: make sure labels exist as cellstr
    if isstring(raw.labels), raw.labels = cellstr(raw.labels); end

end
