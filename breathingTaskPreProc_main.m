
clear

codePre = 'G:\My Drive\GitHub\';
datPre = { ...
    'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ...
    'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\', ...
    'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

% prefix index for data folder:
% 1 = OBE_dupi
% 2 = OBEControl
% 3 = EEG_breathing
datPrei = [ ...
    1,1,1,2,3,3,3,3,3,3,2,1,1,1,1,2,1,1,1,1,1,1,1, ...
    1,1,1,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3];

sessionIDs = { ...
    '250818_Dupi_NMH_JH_1', ...
    '250623_DUPI_NMH_KS_2', ...
    '250623_Dupi_NMH_KS_1', ...
    '250908_OBE_NWU_AS', ...
    '250723_EEG_NWU_IN', ...
    '250725_EEG_NWU_BN', ...
    '250815_EEG_NWU_PP', ...
    '250819_EEG_NWU_ZL', ...
    '250723_EEG_NWU_BK', ...
    '250912_EEG_NWU_JN', ...
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
    '251030_Dupi_NMH_DB_1', ...
    '251110_Dupi_NMH_PC_1', ...
    '251030_Dupi_NMH_DB_2', ...
    '251120_Dupi_NMH_JL_1', ...
    '250818_Dupi_NMH_JH_3', ...
    '251003_EEG_NWU_TI', ...
    '251008_EEG_NWU_JC', ...
    '251008_EEG_NWU_GM', ...
    '251009_EEG_NWU_JM', ...
    '251009_EEG_NWU_SM', ...
    '251027_EEG_NWU_AS', ...
    '251105_EEG_NWU_GL', ...
    '251110_EEG_NWU_GA', ...
    '251111_EEG_NWU_VW', ...
    '251113_EEG_NWU_GH', ...
    '251118_EEG_NWU_ADtest', ...
    '251202_EEG_NWU_GJ', ...
    '260109_EEG_NWU_AA', ...
    '251205_EEG_NWU_AK', ...
    '251208_EEG_NWU_ZA'};
                % '250623_Dupi_NMH_KS_3' not included due to photodiode
                % missing issue. Could fix later. 

% %there are multiple respiration channels in many recordings
% %which one is right for each session: 
% rspIDX = [1,1,3,1,1,1,1,1,1,1]; 
% rspFlip = [1,1,-1,1,1,1,1,-1,-1,1]; %hard code flip


addpath(genpath([codePre 'ZelanoLabScripts']))
addpath(genpath([codePre 'slowBreathing']))
addpath(genpath([codePre 'breathmetrics']))
addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))

targTraceDir = 'G:\My Drive\cZelano\breathingDataFiles'; 
figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';

EEGLOC = readtable(fullfile(codePre, 'ZelanoLabScripts','eegLocs_standard_coords.csv'));
set(0, 'defaultfigurewindowstyle', 'normal')
success = ones(length(sessionIDs),1); 
for s = 1:numel(sessionIDs)
    try
    disp(['working on ', sessionIDs{s}])
    % --- Session descriptor (adjust to your system) ---
    S = struct; 
    S.id   = sessionIDs{s};
    S.root = datPre{datPrei(s)};     % holds exampCueTaskDat.mat
    S.fig  = fullfile(figPath, S.id);
    
    % --- Check if the session is already done --- %
    preDir = fullfile(S.root, S.id, 'preProc');

    outDat = load(fullfile(preDir, [S.id '_breathingPreproc.mat']));
    try
        outDat = outDat.out; 
    catch
        outDat = outDat.chanDat; 
    end
    if isfield(outDat, 'baseEmotion')
        disp(['Done with ' S.id ' ; ' num2str(s)])
        continue
    end

    % --- Params + raw load ---
    [raw, P] = getSessionParams_breathingTask(S);
    %% quick debug stopper: 
%     catch
%         disp('oops')
%     end
% end

%% end debug stopper


    disp(['........................Loaded ', sessionIDs{s}])
    % --- Assemble, preprocess shared pieces ---
    outDat = assemble_outDat_breathing_cue_Task(raw, S, P);
    
    outDat = downsample_data(outDat, P.fs_target);

    if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end

    if P.hasMacros, outDat = preprocess_macros(outDat, P); end
    
    
    disp(['........................spike and blink ', sessionIDs{s}])

    outDat = process_respiration_breathing(outDat, P); 
    if isfield(outDat, 'TTL')
        TaskBreaks = [outDat.TTL/outDat.fs size(outDat.data,2)/outDat.fs];
    else
        TaskBreaks = 0:300:max(outDat.behDat.order)*300-10; 
    end
    for cndi = 1:length(TaskBreaks)
       
        outDat.bmObj(outDat.bmObj(:,2)>TaskBreaks(cndi),12) = cndi; 
    
    end
    outDat.moreThan1 = 1; 
    outDat.rspIDX = P.rspIDX;
    outDat.rspFlip = P.rspFlip; 

    if strcmp(S.id, '250811_Dupi_NMH_TPB_1')
        disp(['TPB 1 had target trace problems'])
    else
        % This doesn't work for the wave breathing task really at all: 
        [outDat, targTraces] = alignTargetBreathingTraceSimplify(outDat, targTraceDir);
    end
    outDat = build_behavior_table_breathingTask(outDat, outDat.bmObj);

    outDat = processECG(outDat, P);
    
    disp(['........................breath behave heart ', sessionIDs{s}])
 
    outDat = flagBadBreaths(outDat); 

    %processing respiration just for plotting
    R = preprocess_respiration_wholetrace(outDat);
    plot_sniff_epochs(outDat, R);

    plotBreathLengths(outDat, R)
   
    writetable(outDat.behDat, [codePre 'closed-loop-respiration\processedBehavior\' ...
                    outDat.sessID  '_processedBreathing.csv']);
   
    % --- Save ---
    preDir = fullfile(S.root, S.id, 'preProc');
    if ~exist(preDir,'dir'), mkdir(preDir); end
    parSave(fullfile(preDir, [S.id '_breathingPreproc.mat']), outDat);
  
    catch ME
        success(s) = 0; 
        disp(['fail for ', sessionIDs{s}, ': ', ME.message])
    end
    
    
end

% USE THIS CODE TO REMOVE AND REDO MACRO CLEANING TO ADJUST SPIKE THRESHOLD
% % labels we want to remove
% rmLabs = {'macBP1','macBP2','macBP3','macBP4','macBP5','spikeCleanVec'};
% 
% % Convert labels to string array for easy comparison
% lblStr = string(outDat.labels);
% 
% % Indices of rows to remove
% rmIdx = ismember(lblStr, rmLabs);
% 
% % Remove rows from data and labels
% outDat.data(rmIdx, :)   = [];
% outDat.labels(rmIdx)    = [];
% 
% % Reset spikeRemoval flag
% outDat.spikeRemoval = 0;
