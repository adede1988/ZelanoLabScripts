
clear
codePre = 'C:\Users\Adam\Documents\GitHub\';

addpath(genpath('C:\Users\Adam\Documents\eeglab2026.0.0'))
addpath(genpath([codePre 'ZelanoLabScripts']))

figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';
% Load once, reuse
EEGLOC = readtable([codePre 'ZelanoLabScripts/eegLocs_standard_coords.csv']);

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
  % <<< all subject-specifics here
  P             = applyParams('O15', S.id);
  disp(['........................Loaded ', sessionIDs{s}])
  % trialStarts, buttonPresses, sniffMarks
  [outDat, raw, TTL] = assemble_outDat_all(S, P);   % loads raw, detects TTLs, assembles

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

  


  outDat.behDat = build_behavior_table_O15(sniffs, raw.beh);

  outDat = refine_onsets_with_phase(outDat, R, P); % uses precomputed phase
  
  plot_sniff_epochs(outDat, R);
  disp(['........................breath behave ', sessionIDs{s}])
  if ~isfolder(fullfile(outDat.OGdataDir,  'preProc'))
      mkdir(fullfile(outDat.OGdataDir,  'preProc'))
  end
  save(fullfile(outDat.OGdataDir,  'preProc', ...
                [outDat.sessID '_O15preproc.mat']), ...
                'outDat', "-v7.3")
  else
       disp(['finished file detected for: ', sessionIDs{s}])
  end
end

