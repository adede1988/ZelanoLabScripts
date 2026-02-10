
clear

codePre = 'G:\My Drive\GitHub\';
datPre = { 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\',...
    'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl'};

outPre = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror';

allPreProcFiles = []; 
%% find preproc files: 
for ii = 1:length(datPre)

    subdir = dir(datPre{ii}); 
    splitNames = cellfun(@(x) strsplit(x, '_'), {subdir.name}, ...
        'uniformoutput', false);

    nameLen = cellfun(@(x) length(x)>3, splitNames);
    isNum = cellfun(@(x) ~isempty(str2num(x{1})), splitNames); 

    subdir = subdir(nameLen & isNum);

    for jj = 1:length(subdir)
        innerDir = dir([subdir(jj).folder '/' subdir(jj).name]); 
        innerDir = innerDir(cellfun(@(x) strcmp('preProc', x), ...
                        {innerDir.name}));
        if length(innerDir) == 1
            preProcDir = dir([innerDir.folder '/' innerDir.name]); 
            idx = cellfun(@(x) ~contains(x, 'macPowPhase'), ...
                {preProcDir.name});
            preProcDir = preProcDir(idx);
            preProcDir = preProcDir(cellfun(@(x) contains(x, '.mat'), ...
                {preProcDir.name}));
            if isempty(allPreProcFiles)
                allPreProcFiles = preProcDir;
            else
                allPreProcFiles = vertcat(allPreProcFiles, preProcDir); 
            end
        else
            warning(['wrong number of preProc folders: ' subdir(jj).name])
        end
    end



end


%% just going after breathing task right now

test = cellfun(@(x) ~contains(x, 'breathing'), {allPreProcFiles.name});
allPreProcFiles(test,:) = [];


%%

addpath([codePre 'ZelanoLabScripts'])
set(0, 'defaultfigurewindowstyle', 'docked')



for ii = 3:length(allPreProcFiles)
    try
    outDat = load([allPreProcFiles(ii).folder '/'...
                     allPreProcFiles(ii).name]);
    
    outDat = outDat.out; 
    
    %ALL TASKS MUST HAVE: 
    chanTmp = struct; 
    %task 
    chanTmp.task = outDat.task; 
    if strcmp(chanTmp.task, 'breathing')
        chanTmp.task = 'breathingTask';
    end
    %subID, sessID, sessNum
    chanTmp.sessID = outDat.sessID; 
    if strcmp(chanTmp.sessID, '250811_Dupi_NMH_TPB_1')
        chanTmp.sessID = '250811_Dupi_NMH_TB_1';
    end
    nameBits = strsplit(chanTmp.sessID, '_'); 
    chanTmp.subID = nameBits{4};
    if length(nameBits) == 5
        chanTmp.sessNum = str2num(nameBits{5}); 
    end
    %OGdataDir
    chanTmp.OGdataDir = strsplit(allPreProcFiles(ii).folder, 'preProc'); 
    chanTmp.OGdataDir = chanTmp.OGdataDir{1}; 
    %loadFile
    tmp = dir(chanTmp.OGdataDir);
    idx = cellfun(@(x) contains(x, 'LoadData'), {tmp.name});
    if sum(idx) == 1
        chanTmp.loadFile = tmp(idx).name; 
    else
        try
            ADidx = cellfun(@(x) contains(x, 'AD'), {tmp(idx).name});
            idx = find(idx); 
            chanTmp.loadFile = tmp(idx(ADidx)).name;
        catch
            error('unidentified load file')
        end
    end

    %type (EEG, OBE, Dupi)
    chanTmp.type = outDat.type; 
    %age
    %sex
    %fs
    chanTmp.fs = outDat.fs; 
    
    %elecLabels (human readable labels)
    chanTmp.labels = outDat.labels;
    %spikeRemoval (1 = yes; 0 = no)
    if isfield(outDat, 'spikeRemoval')
        chanTmp.spikeRemoval = outDat.spikeRemoval; 
    else
        chanTmp.spikeRemoval = 0; 
        warning(['no spike removal: ' chanTmp.subID])
    end
    %eyeBlinkRemoval (1 = yes; 0 = no)
    if isfield(outDat, 'blinkRemoval')
        chanTmp.eyeBlinkRemoval = outDat.blinkRemoval; 
    else
        chanTmp.eyeBlinkRemoval = 0; 
    end
    %spikeEyeFileDir (where is the 2 x time x trial info on 
    %                   spike and blink detections)
    %spikeEyeFileName 
    idx1 = cellfun(@(x) strcmp(x, 'spikeCleanVec'), outDat.labels);
    idx2 = cellfun(@(x) strcmp(x, 'blinkIndicator'), outDat.labels);
    idx3 = cellfun(@(x) strcmp(x, 'badTS'), outDat.labels);
    if sum(idx1+idx2) > 0
        chanTmp.QCFileDir = [outPre '/cleanFiles']; 
        chanTmp.QCFileName = [chanTmp.sessID '_' 'cleaningVecs.mat']; 
        cleanDat = chanTmp; 
        tmp = outDat.data(idx1==1,:); 
        cleanDat.spikeCleanVec = tmp; 
        tmp = outDat.data(idx2==1,:); 
        cleanDat.blinkCleanVec = tmp; 
        tmp = outDat.data(idx3==1,:); 
        cleanDat.badTS = tmp; 
        save([chanTmp.QCFileDir '/' chanTmp.QCFileName],...
                                        'cleanDat'); 
        
    else
        disp([chanTmp.subID ' has no QC information'])
        chanTmp.QCFileDir = []; 
        chanTmp.QCFileName = []; 
    end
    %behDat (task specific behavioral data; grab trial X column matrix
    %task specific
    switch chanTmp.task
        case 'O15'
            chanTmp.behDat = outDat.behDat; 
        case 'PEA_threshold'
            chanTmp.behDat = outDat.behDat; 
        case 'breathingTask'
            chanTmp.behDat = outDat.behDat; 
            chanTmp.baseEmotion = outDat.baseEmotion;
        case 'cueTask'
            chanTmp.behDat = outDat.behDat; 
        otherwise
            error('unknown behavior')
    end  
   

    %respiration
    
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx,:); 
    rspDat = rspDat(outDat.rspIDX,:);
    rspDat = rspDat .* outDat.rspFlip;

    chanTmp.rsp = rspDat; 
    


    %task specific extra stuff: 
    switch chanTmp.task
        case 'O15'
            %rawBehavior: trials X 3 matrix of behavioral responses
            error('check extra stuff for O15')
        case 'PEA_threshold'
            error('check extra stuff for PEA thresh')
        case 'breathingTask'
            idx = cellfun(@(x) contains(x, 'RRint'), outDat.labels);
            chanTmp.RRint = squeeze(outDat.data(idx, :)); 
            idx = cellfun(@(x) contains(x, 'targTrace'), outDat.labels);
            chanTmp.targTrace = squeeze(outDat.data(idx, :)); 
        case 'cueTask'
            error('check extra stuff for O15')
        otherwise
            error('unknown extra')
    end  

   eegChans = {'Fp1', 'Fz' , 'F3' , 'F7' , 'FT9', 'FC5', 'FC1', ...
                'C3' , 'T7' , 'TP9', 'CP5', 'CP1', 'Pz' , 'P3' , ...
                'P7' , 'O1' , 'Oz' , 'O2' , 'P4' , 'P8' , 'TP10',...
                'CP6', 'CP2', 'Cz' , 'C4' , 'T8' ,'FT10', 'FC6', ...
                'FC2', 'F4' , 'F8' , 'Fp2'};
    
    for chi = 1:length(outDat.labels)
        chanDat = chanTmp; 
        lab = chanDat.labels{chi}; 
        if contains(lab, 'macBP') || ...
           sum(cellfun(@(x) strcmp(x, lab), eegChans))==1 
            %chi (index of current channel)
            chanDat.chi = chi; 
            %trialDat (time X trials)
            chanDat.data = squeeze(outDat.data(chi,:)); 
            %chanType (EEG, macro)
            if contains(lab, 'mac')
                chanDat.chanType = 'macro'; 
            else
                chanDat.chanType = 'EEG';
                
            end

            save([outPre '/CHANDAT/' ...
                chanDat.sessID '_' chanDat.chanType '_' chanDat.task '_' ...
                num2str(chi) '.mat'], 'chanDat');
        end
    end
    catch
        disp(['error on: ' allPreProcFiles(ii).name])
    end
end



