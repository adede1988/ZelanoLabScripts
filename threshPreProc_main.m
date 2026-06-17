
clear

codePre = 'C:\Users\Adam\Documents\GitHub\';

addpath(genpath([codePre 'ZelanoLabScripts']))
addpath(genpath('C:\Users\Adam\Documents\eeglab2026.0.0'))


figPath = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\';

EEGLOC = readtable(fullfile(codePre, 'ZelanoLabScripts','eegLocs_standard_coords.csv'));

set(0, 'defaultfigurewindowstyle', 'normal')

cfg        = applyParams('threshTask','main');
sessionIDs = cfg.sessionIDs;

for s = 1:numel(sessionIDs)
    % --- Session descriptor (adjust to your system) ---
    S.id   = sessionIDs{s};
    S.root = cfg.root{s};     % holds exampCueTaskDat.mat
    S.fig  = fullfile(figPath, S.id);
    disp(['working on ', sessionIDs{s}])
    preDir = fullfile(S.root, S.id, 'preProc');
    outDat = load(fullfile(preDir, [S.id '_PEA_threshold_preproc.mat']));
    
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
    P = applyParams('threshTask', S.id);

    % --- Assemble, preprocess shared pieces ---
    [outDat, raw, TTL] = assemble_outDat_all(S, P); %this works for thresh task too!
     disp(['........................Loaded ', sessionIDs{s}])
  % trialStarts, buttonPresses, sniffMarks    
    outDat = downsample_data(outDat, P.fs_target);
    if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end
    outDat = preprocess_macros(outDat, P);
    
  disp(['........................spike and blink ', sessionIDs{s}])
    R = preprocess_respiration_wholetrace(outDat); % fields: rsp, rsp_smooth, phase, onset_metric

    
    trial = 1:45; 
    sniff = outDat.TTL.sniff; 
    x = table(sniff(:)-1000, trial(:), sniff(:), ...
        'variablenames', {'start', 'trial', 'sniff'}); 
    outDat.TTL = x; 
    sniffs = detect_sniffs_from_TTLs(R, P, outDat);  % returns table or matrix
    
    %there's more than one sniff per trial
    
    outDat.rspIDX = P.rspIDX;
    outDat.rspFlip = P.rspFlip; 


    outDat.behDat = build_behavior_table_threshTask(sniffs, raw.beh);

    outDat = refine_onsets_with_phase(outDat, R, P); % uses precomputed phase


    plot_sniff_epochs(outDat, R);
   outDat.moreThan1 = 0; 
  disp(['........................breath behave ', sessionIDs{s}])
    % --- Save ---
    preDir = fullfile(S.root, S.id, 'preProc');
    if ~exist(preDir,'dir'), mkdir(preDir); end
    save(fullfile(preDir, [S.id '_PEA_threshold_preproc.mat']), 'outDat','-v7.3');
    
end

