
clear
% ---- machine paths (everything machine-specific comes from labPaths) ----
addpath(fileparts(mfilename('fullpath')));   % ensure labPaths() is reachable
L       = labPaths();
codePre = L.codePre;
addpath(genpath(L.repo))
addpath(genpath(L.eeglab))

figPath = L.figPath;
EEGLOC  = readtable(L.eegLocCsv);   % load once, reuse

% =====================================================================
%  O15 preprocessing -- main pipeline
%  Sections marked TASK-SHARED are identical across all four pipelines
%  (breathing / cue / thresh / O15); do NOT edit them when adding a new task.
%  Sections marked TASK-SPECIFIC must be rewritten per task.
%  TASK-SPECIFIC pieces for O15 (rewrite these for a new task):
%    - assemble_outDat_all.m  (O15 branch: raw_O15 load + detect_ttls_O15)
%    - detect_ttls_O15.m      (photodiode -> TTL table)
%    - build_behavior_table_O15.m
%  Everything else is SHARED: applyParams, downsample_data, preprocess_eeg,
%  preprocess_macros, preprocess_respiration_wholetrace, detect_sniffs_from_TTLs,
%  refine_onsets_with_phase, paramCheck, writeParams, writePreProcX, plot_sniff_epochs.
% =====================================================================

cfg        = applyParams('O15','main');
sessionIDs = cfg.sessionIDs;

for s = 1:numel(sessionIDs)

  

  disp(['working on ', sessionIDs{s}])
  S.id   = sessionIDs{s};
  S.root = cfg.root{s};
   matPath = fullfile(S.root, S.id);
  if ~exist(fullfile(matPath,  'preProc', ...
                [S.id '_O15preproc.mat']), 'file')
  S.figPath = figPath; 
  if ~isfolder(fullfile(S.figPath, S.id))
      mkdir(fullfile(S.figPath, S.id))
  end
  P             = applyParams('O15', S.id);


  disp(['........................Loaded ', sessionIDs{s}])
  % trialStarts, buttonPresses, sniffMarks
  [outDat, raw, TTL] = assemble_outDat_all(S, P);   % loads raw, detects TTLs, assembles

  if strcmp(P.paramSource, 'guess')
        [outDat, P] = paramCheck(outDat, P);
  end

  outDat.TTL = TTL;

  
  outDat = downsample_data(outDat, P.fs_target);

  if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end


  outDat = preprocess_macros(outDat, P);
  disp(['........................spike and blink ', sessionIDs{s}])
    
  %there's more than one sniff per trial
  outDat.moreThan1 = 1; 
  outDat.rspIDX = P.rspIDX;
  outDat.rspFlip = P.rspFlip; 

  % Precompute whole-trace respiration features ONCE
  R = preprocess_respiration_wholetrace(outDat); % fields: rsp, rsp_smooth, phase, onset_metric

  sniffs = detect_sniffs_from_TTLs(R, P, outDat);  % returns table or matrix

    if strcmp(P.paramSource, 'guess')
        error('check that onsets have been well-detected')
        
    end
    P.paramSource = 'curated'; 
    writeParams(P, S.id);


  % ----- TASK-SPECIFIC (O15): behavior table from sniffs + raw behavior -----
  outDat.behDat = build_behavior_table_O15(sniffs, raw.beh);
  % ----- end TASK-SPECIFIC -----

  outDat = refine_onsets_with_phase(outDat, R, P); % uses precomputed phase
  
  plot_sniff_epochs(outDat, R);
  disp(['........................breath behave ', sessionIDs{s}])
  if ~isfolder(fullfile(outDat.OGdataDir,  'preProc'))
      mkdir(fullfile(outDat.OGdataDir,  'preProc'))
  end
  save(fullfile(outDat.OGdataDir,  'preProc', ...
                [outDat.sessID '_O15preproc.mat']), ...
                'outDat', "-v7.3")

  writePreProcX(P, S.id)   % mark Data Preprocessed = X in dataTracking.xlsx
  else
       disp(['finished file detected for: ', sessionIDs{s}])
  end

end

