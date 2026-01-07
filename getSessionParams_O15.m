function [raw, P] = getSessionParams_O15(S)
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

    matPath = fullfile(root, sessID, 'raw', 'raw_O15', 'raw_O15.mat');
    if ~exist(matPath, 'file')
        error('load_session:MissingMat', 'Raw MAT not found: %s', matPath);
    end
    dat = load(matPath);
    if isfield(dat, 'curDat')
        dat = dat.curDat;
    else
        error('load_session:BadMat', 'Expected curDat in MAT file: %s', matPath);
    end

    behPath = fullfile(root, sessID, 'Behavioral_data', 'O15', ...
                       sprintf('O15_responses_%s.csv', sessID));
    if ~exist(behPath, 'file')
        error('load_session:MissingBehavior', 'Behavior CSV not found: %s', behPath);
    end
    behTbl = readtable(behPath);

   

    % Populate raw struct
    raw = struct();
    raw.sessID    = sessID;
    raw.fs_raw    = dat.rawData.fsample;
    raw.data      = dat.rawData.trial{1};
    raw.labels    = dat.outLabs;
    if isfield(dat, 'ncslabels'), raw.ncslabels = dat.ncslabels; end
    raw.beh       = behTbl;
    raw.paths.mat = matPath;
    raw.paths.beh = behPath;
    raw.paths.fig = fullfile(S.figPath, S.id);

    % Heuristic type (kept for bookkeeping)
    if contains(sessID, 'OBE', 'IgnoreCase', true) || contains(root, 'OBE', 'IgnoreCase', true)
        raw.type = 'OBE';
    else
        raw.type = 'Dupi';
    end

    %% --- Defaults (session-agnostic) ---
    P = struct();
    P.task      = 'O15';
    P.type      = raw.type;
    P.fs_target = 500;          % target downsample Fs
    P.debug     = false;        % turn plots on/off upstream

    % Resp channel selection (from your arrays: all 1's, all no-flip)
    P.rspIDX    = 1;
    P.rspFlip   = 1;

    % Respiration onset metric thresholds (global defaults)
    P.respThresh   = 500;       % default 'test' threshold
    P.cuedBackBuff = 150;       % samples (downsampled space)
    P.adjWin       = 500;       % samples for backtrack during refine

    % EEG / spike cleaning defaults
    P.hasEEG      = false;
    P.spikeThresh = 20;
    P.spikeWin    = 11;
    P.spikeClean  = false;
    P.macroRemove = []; 

    % Photodiode pulse parsing 
    P.pd = struct();
    P.pd.zthresh        = -2;    % z-score threshold for "low"
    P.pd.minPulseSamp   = 200;   % remove pulses < 200 samples
    P.pd.maxPulseSamp   = 1200;  % remove pulses > 1200 samples
    P.pd.trialSplitSamp = 850;   % len < 850 -> trial marks; > 850 -> sniffs

    % TTL expectations / overrides
    P.ttl = struct();
    P.ttl.expectedTrialCount   = 30;  % 15 trials × 2 marks (start & button)
    P.ttl.removeTrialMarksIdx  = [];  % indices to drop from trialMarks after detection (if any)
    P.ttl.note                 = "";  % free-form notes

    %% --- Session-specific overrides (aggregated from your script) ---
    switch sessID
        % EEG present + spike cleaning parameters
        case '250904_OBE_NWU_TI'
            P.hasEEG      = true;
            P.spikeThresh = 50;
            P.spikeWin    = 11;
            P.spikeClean  = true;

            % Resp thresholds
            P.respThresh   = 5000;
            P.cuedBackBuff = 100;

            % Refine window narrower for this participant
            P.adjWin       = 300;

            % Known TTL anomaly in your script
            P.ttl.removeTrialMarksIdx = 28;
            P.ttl.note = "Drop aberrant extra TTL at index 28.";

        case '250623_DUPI_NMH_KS_2'
            P.hasEEG      = true;
            P.spikeThresh = 15;
            P.spikeWin    = 11;
            P.spikeClean  = false;

            P.respThresh   = 3000;
            P.cuedBackBuff = 350;

        case '250623_Dupi_NMH_KS_1'
            P.hasEEG      = false;
            P.spikeThresh = 15;
            P.spikeWin    = 11;
            P.spikeClean  = false;
            P.macroRemove = 6; 
            P.respThresh   = 5000;
            P.cuedBackBuff = 350;
        case '250929_Dupi_NMH_GH_1'
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.macroRemove = 6; 
            P.rspFlip     = 1;

        case '250818_Dupi_NMH_JH_2'
            P.hasEEG      = true;
            P.spikeThresh = 20;
            P.spikeWin    = 11;
            P.spikeClean  = true;

            % Uses defaults for resp thresholds

        case '250818_Dupi_NMH_JH_1'
            P.hasEEG      = true;
            P.spikeThresh = 20;
            P.spikeWin    = 11;
            P.spikeClean  = true;

        case '250908_OBE_NWU_AS'
            P.hasEEG      = true;
            P.spikeThresh = 20;
            P.spikeWin    = 11;
            P.spikeClean  = false;

            P.respThresh   = 3000;
            P.cuedBackBuff = 150;

        case '250811_Dupi_NMH_TPB_1'
            P.hasEEG      = true;
            P.spikeThresh = 10;
            P.spikeWin    = 9;
            P.spikeClean  = true;

            P.respThresh   = 20;
            P.cuedBackBuff = 150;
        case '250811_Dupi_NMH_TB_2'
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.rspFlip     = -1;
                
            P.respThresh   = 20;
            P.cuedBackBuff = 150;
        case '251002_Dupi_NMH_AB_1'
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.respThresh   = 2000;
        case '251006_OBE_NWU_RY_1'
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.respThresh   = 4000;
            P.cuedBackBuff = 400;
        case '251027_Dupi_NMH_DL_1'
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.cuedBackBuff = 300;
            P.respThresh   = 4000;
        case '251013_Dupi_NMH_JN_2'
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.macroRemove = 6;
            P.cuedBackBuff = 200;
            P.respThresh   = 2000;
        case '251009_OBE_NWU_CP_1' 
            P.macroRemove = 6; 
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.respThresh   = 2000;
            P.cuedBackBuff = 250;
        case '250929_Dupi_NMH_GH_2'
            P.rspIDX      = 1; 
            P.macroRemove = [5,6]; 
            P.hasEEG      = true;
            P.spikeClean  = true;
            P.respThresh   = 4000;
            P.cuedBackBuff = 350;
            P.ttl.removeTrialMarksIdx = [27,28];
            P.ttl.note = "Drop extra TTLs from restarting task.";
        case '251002_Dupi_NMH_AB_2'
            P.rspIDX      = 1; 
            P.macroRemove = [4,5,6]; 
            P.hasEEG      = true;
            P.spikeClean  = false;
            P.respThresh   = 1000;
        otherwise
            % keep defaults
    end

    % sanity: make sure labels exist as cellstr
    if isstring(raw.labels), raw.labels = cellstr(raw.labels); end

end
