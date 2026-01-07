
clear

codePre = 'G:\My Drive\GitHub\';
datPre = { 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ... 
           'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\'};

%prefix index for data folder: 
datPrei = [1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,1]; 

sessionIDs = {'250818_Dupi_NMH_JH_1', ... %preprocessed
               '250623_DUPI_NMH_KS_2',... 
               '250623_Dupi_NMH_KS_1',... 
                '250818_Dupi_NMH_JH_2',... 
                '250811_Dupi_NMH_TPB_1',... 
              '230611_OBE_NMH_AZ',... 
                '241017_OBE_NMH_AS',...  
                '240923_OBE_NMH_HRM',... 
                '250310_OBE_NMH_FS',...
                '250313_OBE_NMH_CS',...
                 '250929_Dupi_NMH_GH_1',...
            '251002_Dupi_NMH_AB_1',...
            '251027_Dupi_NMH_DL_1', ...
            '250929_Dupi_NMH_GH_2',...
            '251002_Dupi_NMH_AB_2'...
            '251013_Dupi_NMH_JN_2'};  

% %there are multiple respiration channels in many recordings
% %which one is right for each session: 
% rspIDX = [1,1,3,1,1,1,1,1,1,1]; 
% rspFlip = [1,1,-1,1,1,1,1,-1,-1,1]; %hard code flip


addpath(genpath([codePre 'ZelanoLabScripts']))
addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))


figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';

EEGLOC = readtable(fullfile(codePre, 'ZelanoLabScripts','myEEGcoords_thetaPhi.csv'));


for s = 1:numel(sessionIDs)
    % --- Session descriptor (adjust to your system) ---
    S.id   = sessionIDs{s};
    S.root = datPre{datPrei(s)};     % holds exampCueTaskDat.mat
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
    [raw, P] = getSessionParams_cueTask(S);
    
    % --- Assemble, preprocess shared pieces ---
    outDat = assemble_outDat_breathing_cue_Task(raw, S, P);
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
    outDat.rspIDX = P.rspIDX;
    outDat.rspFlip = P.rspFlip; 


    outDat.behDat = build_behavior_table_cueTask(sniffs, raw.beh);

    outDat = refine_onsets_with_phase(outDat, R, P); % uses precomputed phase


    plot_sniff_epochs(outDat, R);
   
  disp(['........................breath behave ', sessionIDs{s}])
    % --- Save ---
    preDir = fullfile(S.root, S.id, 'preProc');
    if ~exist(preDir,'dir'), mkdir(preDir); end
    save(fullfile(preDir, [S.id '_cueTaskPreproc.mat']), 'outDat','-v7.3');
    
end

