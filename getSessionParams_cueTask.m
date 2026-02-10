function [raw, P] = getSessionParams_cueTask(S)
% getSessionParams_cueTask
%   Loads the prebuilt outDat for the cue task and aggregates
%   all subject-specific parameters into P (centralizing what
%   was hard-coded in cueTaskPreProc_scratch.m).
%
% Inputs
%   S.id   : session ID (e.g., '250904_OBE_NWU_TI')
%   S.root : folder containing <S.id>.mat (with outDat inside)
%
% Outputs
%   raw : struct similar to O15 loader outputs
%         .sessID, .fs_raw, .data, .labels, .beh, .TTL (if present)
%   P   : parameters for preprocessing and event handling
%         .task, .type, .fs_target, .debug
%         .computeResp, .rspIDX, .rspFlip
%         .hasEEG, .spikeThresh, .spikeWin, .spikeClean
%         .respThresh, .cuedBackBuff, .adjWin
%         .ttlMap (names to look for if TTL is a table/struct)

    arguments
        S struct
    end

    %% --- Load prebuilt outDat ---
    matPath = fullfile(S.root, S.id, 'preProc', [S.id '_cueTaskPreProc.mat']);
    if ~exist(matPath, 'file')
        error('getSessionParams_cueTask:MissingMat', ...
              'Expected MAT at %s (containing outDat).', matPath);
    end
    tmp = load(matPath);
   
    if ~isfield(tmp, 'outDat')
        error('getSessionParams_cueTask:NoOutDat', ...
              'MAT must contain variable outDat.');
    end
    od = tmp.outDat;

    % raw view to match O15 pipeline expectations
    raw = struct();
    raw.sessID = char(od.sessID);
    raw.fs_raw = od.fs;
    raw.data   = od.data;
    raw.labels = od.labels;
    raw.beh    = od.behDat;
    if isfield(od, 'TTL'), raw.TTL = od.TTL; end
    raw.paths  = struct('root', S.root);

    sessID = char(S.id);
    if contains(sessID, 'OBE', 'IgnoreCase', true) || contains(char(S.root), 'OBE', 'IgnoreCase', true)
        raw.type = 'OBE';
    else
        raw.type = 'Dupi';
    end
    %% --- Defaults (session-agnostic) ---
    P = struct();
    P.task        = 'cueTask';
    P.type        = raw.type;
    P.fs_target   = 500;
    P.debug       = false;

    % Respiration (optional — leave on so downstream helpers work)
    P.computeResp = true;
    P.rspIDX      = 1;    
    P.rspFlip     = 1;    

    % EEG / spike cleaning defaults
    P.hasEEG      = true;
    P.spikeThresh = 20;
    P.spikeWin    = 11;
    P.spikeClean  = true;
    P.macroRemove = []; 

    % Respiration onset metric & alignment defaults
    P.respThresh   = 500;   % scratch "otherwise" default
    P.cuedBackBuff = 150;   % scratch "otherwise" default
    P.adjWin       = 500;   % scratch default, except TI case

    % TTL name mapping (used if raw.TTL is a struct/table)
    P.ttlMap = struct( ...
        'cue',    {'cue','Cue','cueOnset'}, ...
        'target', {'targ','target','TargetOnset'}, ...
        'resp',   {'resp','response','button'} );

    %% --- Subject-specific overrides from cueTaskPreProc_scratch.m ---

     switch raw.sessID

        case '250818_Dupi_NMH_JH_1'
            P.rspIDX = 1;  P.rspFlip = 1;
            P.hasEEG = true;  P.spikeClean = true;
            P.respThresh  = 4000; P.cuedBackBuff = 350;

        case '250623_DUPI_NMH_KS_2'
            P.rspIDX = 1;  P.rspFlip = 1;
            P.hasEEG = true;  P.spikeClean = false;
            P.respThresh  = 3000; P.cuedBackBuff = 350;
        
        case '250623_Dupi_NMH_KS_1'
            P.rspIDX = 3;  P.rspFlip = -1;
            P.hasEEG = false; P.spikeClean = false;
            P.macroRemove = 6;
            P.respThresh  = 5000; P.cuedBackBuff = 350;
        
        case '250818_Dupi_NMH_JH_2'
            P.rspIDX = 1;  P.rspFlip = 1;
            P.hasEEG = true;  P.spikeClean = true;
            P.respThresh  = 4000; P.cuedBackBuff = 350;
        
        case '250811_Dupi_NMH_TPB_1'
            P.rspIDX = 1;  P.rspFlip = 1;
            P.hasEEG = true;  P.spikeClean = true;
            P.respThresh  = 20;  P.cuedBackBuff = 150;
        
        case '230611_OBE_NMH_AZ'
            P.rspIDX = 1;  P.rspFlip = 1;
            P.hasEEG = false; P.spikeClean = true;
            P.respThresh  = 3000; P.cuedBackBuff = 150;
        
        case '241017_OBE_NMH_AS'
            P.rspIDX = 1;  P.rspFlip = 1;
            P.hasEEG = false; P.spikeClean = false;
            P.respThresh  = 3000; P.cuedBackBuff = 150;
        
        case '240923_OBE_NMH_HRM'
            P.rspIDX = 1;  P.rspFlip = 1;
            P.hasEEG = false; P.spikeClean = true;
            P.respThresh  = 1000; P.cuedBackBuff = 350;
        
        case '250310_OBE_NMH_FS'
            P.rspIDX = 1;  P.rspFlip = -1;
            P.hasEEG = false; P.spikeClean = true;
            P.respThresh  = 3000; P.cuedBackBuff = 150;
        
        case '250313_OBE_NMH_CS'
            P.rspIDX = 1;  P.rspFlip = -1;
            P.hasEEG = false; P.spikeClean = false;
            P.respThresh  = 3000; P.cuedBackBuff = 150;
        
        case '250929_Dupi_NMH_GH_1'
            P.rspFlip = 1;
            P.macroRemove = 6; P.cuedBackBuff = 250;
            P.hasEEG = true;  P.spikeClean = true;
        
        case '250811_Dupi_NMH_TB_2'
            P.rspFlip = -1;
            P.hasEEG = true;  P.spikeClean = true;
            P.respThresh  = 20;  P.cuedBackBuff = 150;
        
        case '251002_Dupi_NMH_AB_1'
            P.hasEEG = true;  P.spikeClean = true;
            P.respThresh  = 2000; P.cuedBackBuff = 250;
        
        case '251027_Dupi_NMH_DL_1'
            P.hasEEG = true;  P.spikeClean = true;
            P.respThresh  = 4000; P.cuedBackBuff = 300;
        
        case '251013_Dupi_NMH_JN_2'
            P.macroRemove = 6;
            P.hasEEG = true;  P.spikeClean = true;
            P.respThresh  = 2000; P.cuedBackBuff = 200;
        
        case '250929_Dupi_NMH_GH_2'
            P.rspIDX = 1;
            P.macroRemove = [5,6];
            P.hasEEG = true;  P.spikeClean = true;
            P.respThresh  = 4000; P.cuedBackBuff = 350;
        
        case '251002_Dupi_NMH_AB_2'
            P.rspIDX = 1;
            P.macroRemove = [4,5,6];
            P.hasEEG = true;  P.spikeClean = false;
            P.respThresh  = 1000;


        otherwise
            % Keep defaults
    end

end
