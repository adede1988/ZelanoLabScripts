function O15PreProc_main
  codePre = 'G:\My Drive\GitHub\';
  datPre  = {'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ...
             'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\', ...
             'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

  sessionIDs = {'250818_Dupi_NMH_JH_1', ... %preProcessed
               '250623_DUPI_NMH_KS_2',... 
               '250623_Dupi_NMH_KS_1',... 
               '250908_OBE_NWU_AS', ...  
                '250904_OBE_NWU_TI', ...  
                '250818_Dupi_NMH_JH_2',...
                '250811_Dupi_NMH_TPB_1'};           
  datPrei = [1,1,1,2,2,1,1];       

  addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))
  addpath(genpath([codePre 'ZelanoLabScripts']))

  figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';
  % Load once, reuse
  EEGLOC = readtable([codePre 'ZelanoLabScripts/myEEGcoords_thetaPhi.csv']);

  for s = 1:numel(sessionIDs)
      S.id   = sessionIDs{s};
      S.root = datPre{datPrei(s)};
      S.figPath = figPath; 
      % <<< all subject-specifics here
      [raw, P]      = getSessionParams_O15(S);

      % trialStarts, buttonPresses, sniffMarks    
      TTL    = detect_ttls_O15(raw, P);

      outDat = assemble_outDat_O15(raw, S, P);   % copies metadata, sets task/type

      
      outDat = downsample_data(outDat, P.fs_target);

      if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end


      outDat = preprocess_macros(outDat, P);


      % Precompute whole-trace respiration features ONCE
      R = preprocess_respiration_wholetrace(outDat, P); % fields: rsp, rsp_smooth, phase, onset_metric

      sniffs = detect_sniffs_from_TTLs(R, TTL, P, outDat);  % returns table or matrix

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
  end
end
