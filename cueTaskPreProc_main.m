
clear
% ---- machine paths (everything machine-specific comes from labPaths) ----
addpath(fileparts(mfilename('fullpath')));   % ensure labPaths() is reachable
L       = labPaths();
codePre = L.codePre;
addpath(genpath(L.repo))
addpath(genpath(L.eeglab))

figPath = L.figPath;
EEGLOC  = readtable(L.eegLocCsv);

set(0, 'defaultfigurewindowstyle', 'normal')

% =====================================================================
%  cueTask preprocessing -- main pipeline
%  TASK-SHARED sections are identical across all four pipelines; do NOT edit
%  them when adding a new task. TASK-SPECIFIC sections must be rewritten.
%  TASK-SPECIFIC pieces for cueTask (rewrite these for a new task):
%    - assembleRaw_cueTask.m  (raw load -> raw struct)
%    - cueTask_makeOutDat.m   (photodiode/behavior ingestion -> _cueTaskPreProc.mat)
%    - outMat_to_table.m      (cue behavioral .mat -> table; used by makeOutDat)
%    - build_behavior_table_cueTask.m
%  Everything else is SHARED: applyParams, assembleOutDat, downsample_data,
%  preprocess_eeg, preprocess_macros, preprocess_respiration_wholetrace,
%  detect_sniffs_from_TTLs, refine_onsets_with_phase, paramCheck, writeParams,
%  writePreProcX, plot_sniff_epochs.
% =====================================================================

cfg        = applyParams('cueTask','main');
sessionIDs = cfg.sessionIDs;

for s = 1:numel(sessionIDs)
    % --- Session descriptor (adjust to your system) ---
    S.id   = sessionIDs{s};
    S.root = cfg.root{s};
    S.fig  = fullfile(figPath, S.id);
    disp(['working on ', sessionIDs{s}])
    preDir = fullfile(S.root, S.id, 'preProc');
    outDat = load(fullfile(preDir, [S.id '_cueTaskPreproc.mat']));
    
    if isfield(outDat, 'out')
        outDat = outDat.out; 
    elseif isfield(outDat, 'outDat')
        outDat = outDat.outDat; 
    else
        error('unexpected missing field in outDat')
    end
   
    if isfield(outDat, 'moreThan1')
        disp(['Done with ' S.id ' ; ' num2str(s)])
        continue
    end
    % --- Params + raw load ---
    P = applyParams('cueTask', S.id);

    if strcmp(P.paramSource, 'guess')
        [outDat, P] = paramCheck(outDat, P);
    end
    outDat.rspIDX = P.rspIDX;
    outDat.rspFlip = P.rspFlip; 


    % --- Assemble: TASK-SPECIFIC loader + shared assembler ---
    raw    = assembleRaw_cueTask(S);    % <-- TASK-SPECIFIC: edit/replace for a new task
    outDat = assembleOutDat(raw, S, P);  % shared
     disp(['........................Loaded ', sessionIDs{s}])
  % trialStarts, buttonPresses, sniffMarks    
    outDat = downsample_data(outDat, P.fs_target);
    if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end
    outDat = preprocess_macros(outDat, P);
    
  disp(['........................spike and blink ', sessionIDs{s}])
    R = preprocess_respiration_wholetrace(outDat); % fields: rsp, rsp_smooth, phase, onset_metric

    sniffs = detect_sniffs_from_TTLs(R, P, outDat);  % returns table or matrix
    
    %there's more than one sniff per trial
    outDat.moreThan1 = 0; 


    % ----- TASK-SPECIFIC (cueTask): behavior table from sniffs + raw behavior -----
    outDat.behDat = build_behavior_table_cueTask(sniffs, raw.beh);
    % ----- end TASK-SPECIFIC -----

    outDat = refine_onsets_with_phase(outDat, R, P); % uses precomputed phase

    if strcmp(P.paramSource, 'guess')
        error('check that onsets have been well-detected')
        
    end
    P.paramSource = 'curated'; 
    writeParams(P, S.id);

    plot_sniff_epochs(outDat, R);
   
  disp(['........................breath behave ', sessionIDs{s}])
    % --- Save ---
    preDir = fullfile(S.root, S.id, 'preProc');
    if ~exist(preDir,'dir'), mkdir(preDir); end
    save(fullfile(preDir, [S.id '_cueTaskPreproc.mat']), 'outDat','-v7.3');

    writePreProcX(P, S.id);   % mark Data Preprocessed = X in dataTracking.xlsx
end

