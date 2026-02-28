
clear
codePre = 'G:\My Drive\GitHub\';
datPre  = {'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ...
         'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\', ...
         'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

sessionIDs = {'250818_Dupi_NMH_JH_1', ...
              '250623_DUPI_NMH_KS_2', ...
              '250623_Dupi_NMH_KS_1', ...
              '250908_OBE_NWU_AS', ...
              '250904_OBE_NWU_TI', ...
              '250818_Dupi_NMH_JH_2', ...
              '250811_Dupi_NMH_TPB_1', ...
              '250811_Dupi_NMH_TB_2', ...
              '250929_Dupi_NMH_GH_1', ...
              '251009_OBE_NWU_CP_1', ...
              '251002_Dupi_NMH_AB_1', ...
              '251027_Dupi_NMH_DL_1', ...
              '250929_Dupi_NMH_GH_2', ...
              '251002_Dupi_NMH_AB_2', ...
              '251013_Dupi_NMH_JN_2', ...
              '251006_OBE_NWU_RY_1', ...
              '251030_Dupi_NMH_DB_1', ...
              '251110_Dupi_NMH_PC_1', ...
              '250623_Dupi_NMH_KS_3', ...
              '251030_Dupi_NMH_DB_2', ...
              '251120_Dupi_NMH_JL_1', ...
              '250818_Dupi_NMH_JH_3'};
  



datPrei = [1,1,1,2,2,1,1,1,1,2,1,1,1,1,1,2,1,1,1,1,1,1];       

addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))
addpath(genpath([codePre 'ZelanoLabScripts']))

figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';
% Load once, reuse
EEGLOC = readtable([codePre 'ZelanoLabScripts/myEEGcoords_thetaPhi.csv']);

for s = 18:numel(sessionIDs)

  

  disp(['working on ', sessionIDs{s}])
  S.id   = sessionIDs{s};
  S.root = datPre{datPrei(s)};
   matPath = fullfile(S.root, S.id);
  if ~exist(fullfile(matPath,  'preProc', ...
                [S.id '_O15preproc.mat']), 'file')    
  S.figPath = figPath; 
  if ~isfolder(fullfile(S.figPath, S.id))
      mkdir(fullfile(S.figPath, S.id))
  end
  % <<< all subject-specifics here
  [raw, P]      = getSessionParams_O15(S);
  disp(['........................Loaded ', sessionIDs{s}])
  % trialStarts, buttonPresses, sniffMarks    
  [TTL, raw]    = detect_ttls_O15(raw, P);

  outDat = assemble_outDat_O15(raw, S, P);   % copies metadata, sets task/type

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

