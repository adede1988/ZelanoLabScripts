clear

codePre = 'G:\My Drive\GitHub\';
datPre = { 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ... 
           'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\'};
behDatPath = 'R:\Neurology\Zelano_Lab\Lab_Common\OBE_task_backup\tasks\pea_threshold\results';

%prefix index for data folder: 
datPrei = [1,1,1,1,1,1,1,1,1]; 

sessionIDs = {'250818_Dupi_NMH_JH_1', ... 
               '250623_DUPI_NMH_KS_2',... 
               '250623_Dupi_NMH_KS_1',... 
                '250818_Dupi_NMH_JH_2',...
                '250811_Dupi_NMH_TPB_1',...
                '250811_Dupi_NMH_TB_2',...
                '250929_Dupi_NMH_GH_1',...
                '251002_Dupi_NMH_AB_1',...
                '251027_Dupi_NMH_DL_1'...
                };   

%there are multiple respiration channels in many recordings
%which one is right for each session: 
rspIDX = [1,1,3,1,1,1,1,1,1]; 
rspFlip = [1,1,-1,1,1,-1,1,1,1]; %hard code flip

addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
addpath(genpath([codePre 'myFrequentUse']))
addpath([codePre 'myFrequentUse/export_fig_repo'])
addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))

addpath([codePre 'fieldtrip-20230118'])
addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])

set(0, 'defaultfigurewindowstyle', 'docked')
ft_defaults



for sessi = 1:length(sessionIDs)
%% data load 
sessi
    % % %find behavioral folder: 
    % % subFolders = dir([datPre{datPrei(sessi)} sessionIDs{sessi}]);
    % % subFolders = subFolders([subFolders.isdir]);
    % % idx = cellfun(@(x) contains(x, 'ehavior'), {subFolders.name});
    % % idx = find(idx); 
    % % if length(idx) ~= 1
    % %     error('behavioral data folder not uniquely identified!')
    % % end
    % % behFold = subFolders(idx);
    datFolders = dir([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\']);

    idx = cellfun(@(x) contains(x, 'raw_PEAintensityPleasantness'), {datFolders.name});
    idx = find(idx); 

    

    
    if isscalar(idx)
        %process when one run: 

        %load in both files: 
        dat = load([datFolders(idx).folder '/' ...
                     datFolders(idx).name  '/' ...
                     datFolders(idx).name  '.mat']);
        dat = dat.curDat; 
     
        

        %find the behavioral files: 
        behDir = dir([behDatPath '/' sessionIDs{sessi}]);
        idx = find(cellfun(@(x) contains(x, '.mat'), {behDir.name})); 

      
        if ~isscalar(idx)
            error('wrong number of behavioral files detected')
        end

        %load the behavioral files
        behDat = [behDir(idx).folder '/' behDir(idx).name];
       
        behDat = load(behDat).outMat;
        behDat = cat(1, behDat{:});
        %fill in empty values
        f = @(x) {735}; 
        emptyIDX = cellfun(@isempty, behDat(:,4));
        behDat(emptyIDX,4) = cellfun(@(x) f(x), behDat(emptyIDX, 4)); 
        emptyIDX = cellfun(@isempty, behDat(:,5));
        behDat(emptyIDX,5) = cellfun(@(x) f(x), behDat(emptyIDX, 5)); 
        behDat = cell2mat(behDat); 
        behDat = array2table(behDat, 'variablenames', ...
                            {'trialNum', 'Odor', 'RT', ...
                             'pleasantness', 'intensity', 'ITI',...
                             'Trial_StartTime', 'Trial_EndTime'});
        
       


        %start making the oudDat struct
        outDat = struct; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
        outDat.sessID = sessionIDs{sessi}; 
        outDat.behDat = behDat;
        
        %process the TTL pulses for dat: 
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat.rawData.trial{1}(idx, :); 

%%

         % Find missing samples in the photodiode channel
        nanidx = isnan(photoDiode);
        nNan   = sum(nanidx);
        
        if nNan > 4000
            error('too many missing values!');
        end
        
        % Get raw data for this trial (rows = channels, cols = time)
        rawData = dat.rawData.trial{1};
        
        if nNan > 0
            'here'
            % --- 1) Interpolate photodiode (1-D vector) ---
            % Fill internal NaNs by linear interpolation
            photoDiode = fillmissing(photoDiode, 'linear');
            % Any leading/trailing NaNs get replaced by nearest neighbor
            photoDiode = fillmissing(photoDiode, 'nearest');
        
            % --- 2) Interpolate rawData along time (dimension 2) ---
            % Linear interpolation across time for each channel
            rawData = fillmissing(rawData, 'linear', 2);
            % Nearest neighbor for any remaining edge NaNs
            rawData = fillmissing(rawData, 'nearest', 2);
        end

        outDat.data = rawData;

        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
    

        switch sessionIDs{sessi}
            case '250623_DUPI_NMH_KS_2'
                diodeThresh = -2; 
            otherwise
                diodeThresh = -1.5;
        end

        downs = find(photoDiode(1:end-1) > diodeThresh & ...
            photoDiode(2:end)<diodeThresh);
        ups = find(photoDiode(1:end-1) < diodeThresh & ...
            photoDiode(2:end)>diodeThresh);
        difVals = ups - downs; %difVals is length of TTL pulses
        downs(difVals>1500) = []; 
        ups(difVals>1500) = []; 
        difVals(difVals>1500) = []; 

        downs(difVals<300) = []; 
        ups(difVals<300) = []; 
        difVals(difVals<300) = []; 

        starti = 1; 
        for di = 3:length(downs)
            if downs(di) - downs(di-2) < 3500
                starti = di; 
            end
        end
        downs(1:starti) = []; 
        difVals(1:starti) = []; 

        if length(downs) ~= 45
            error('wrong number of TTLs')
        end
        
        
    
        figure
        plot(photoDiode)
        xline(downs, 'color', 'magenta', 'linewidth', 2)
    

        %store all TTLs into one matrix: 
        %col 1: sniffs
       
        TTLs = nan(45, 1); 
        TTLs(:,1) = downs; 
        TTLs =round(TTLs ./ 4);

   
        
    
    else
        %process for all in two runs: 
        error('what are you supposed to do?')


    end
    
    outDat.task = "PEAintensityPleasantness_threshold"; 
    outDat.OGdataDir = [datPre{datPrei(sessi)} sessionIDs{sessi}];
    tmp = dir([datPre{datPrei(sessi)} sessionIDs{sessi}]);
    tmp = tmp(cellfun(@(x) contains(x, '.m'), {tmp.name}));
    tmp = tmp(cellfun(@(x) contains(x, 'LoadData'), {tmp.name}));
    if size(tmp,1) == 1
        outDat.loadFile = tmp.name;
    else 
        tmp = tmp(cellfun(@(x) contains(x, 'AD.m'), {tmp.name}));
        if size(tmp,1) == 1
            outDat.loadFile = tmp.name;
        else 
            error('load file not identified uniquely')
        end
    end
    outDat.preProcScript = 'threshPreProc_makeOutDat.m'; 
    if datPrei(sessi) == 1
        outDat.type = 'Dupi'; 
    elseif datPrei(sessi) == 2
        outDat.type = 'OBE';
    end
    outDat.TTL = table; 
    outDat.TTL.sniff = TTLs;


     if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc'], 'dir')
         mkdir([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc']);
     end
    
     save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                        sessionIDs{sessi} '_PEA_threshold_preproc.mat'], ...
                        'outDat', "-v7.3")
   
end

    