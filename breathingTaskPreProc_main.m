
clear

codePre = 'C:\Users\Adam\Documents\GitHub\';

addpath(genpath([codePre 'ZelanoLabScripts']))
addpath(genpath([codePre 'slowBreathing']))
% breathmetrics toolbox is not used in this pipeline -> not loaded
addpath(genpath('C:\Users\Adam\Documents\eeglab2026.0.0'))

targTraceDir = 'G:\My Drive\cZelano\breathingDataFiles'; 
figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';

EEGLOC = readtable(fullfile(codePre, 'ZelanoLabScripts','eegLocs_standard_coords.csv'));
set(0, 'defaultfigurewindowstyle', 'normal')

cfg        = applyParams('breathingTask','main');
sessionIDs = cfg.sessionIDs;

success = ones(length(sessionIDs),1);
for s = 1:numel(sessionIDs)
    try
    disp(['working on ', sessionIDs{s}])
    % --- Session descriptor (adjust to your system) ---
    S = struct;
    S.id   = sessionIDs{s};
    S.root = cfg.root{s};     % holds exampCueTaskDat.mat
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
    P = applyParams('breathingTask', S.id);
    %% quick debug stopper:
%     catch
%         disp('oops')
%     end
% end

%% end debug stopper


    disp(['........................Loaded ', sessionIDs{s}])
    % --- Assemble, preprocess shared pieces ---
    [outDat, raw, TTL] = assemble_outDat_all(S, P);

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
