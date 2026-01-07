
clear
datPre  = {fullfile(pwd, 'sample_data')};

sessionIDs = {'raw_O15'};           
datPrei = [1,1,1,2,2,1,1];       


addpath(genpath(pwd));
%INSERT ADD PATH FOR eeglab2025.0.0 here: 
% addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))

figPath = fullfile(pwd, 'testFigs');
% Load once, reuse
EEGLOC = readtable(fullfile(pwd,'myEEGcoords_thetaPhi.csv'));

%CHAT GPT: IGNORE looping structure and focus on making a working pipeline:
s=1; 
% for s = 1:numel(sessionIDs)
  S.id   = sessionIDs{s};
  S.root = datPre{datPrei(s)};
  S.figPath = figPath; 
  % <<< all subject-specifics here
  [raw, P]      = getSessionParams_O15(S);

  % trialStarts, buttonPresses, sniffMarks    
  TTL    = detect_ttls_O15(raw, P);

  outDat = assemble_outDat_O15(raw, S, P);   % copies metadata, sets task/type

  outDat.TTL = TTL; 

  
  outDat = downsample_data(outDat, P.fs_target);

  if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end


  outDat = preprocess_macros(outDat, P);


  % Precompute whole-trace respiration features ONCE
  R = preprocess_respiration_wholetrace(outDat, P); % fields: rsp, rsp_smooth, phase, onset_metric

  sniffs = detect_sniffs_from_TTLs(R, P, outDat);  % returns table or matrix

  %there's more than one sniff per trial
  outDat.moreThan1 = 1; 
  outDat.rspIDX = P.rspIDX;
  outDat.rspFlip = P.rspFlip; 


  outDat.behDat = build_behavior_table_O15(sniffs, raw.beh);

  outDat = refine_onsets_with_phase(outDat, R, P); % uses precomputed phase
  
  plot_sniff_epochs(outDat, R);

  save(fullfile(outDat.OGdataDir,  'preProc', ...
                [outDat.sessID '_O15preproc.mat']), ...
                'outDat', "-v7.3")
% end

