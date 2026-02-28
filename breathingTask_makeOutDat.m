


%custom function to stitch together breathing recordings into one file

clear

codePre = 'G:\My Drive\GitHub\';
datPre = { 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ... 
           'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\',...
           'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

%prefix index for data folder: 
datPrei = [1,1,1,2,3,3,3,3,3,3,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1]; 

sessionIDs = {'250818_Dupi_NMH_JH_1', ... #preprocessed
               '250623_DUPI_NMH_KS_2',... #preprocessed
               '250623_Dupi_NMH_KS_1',... #preprocessed
               '250908_OBE_NWU_AS', ...
                '250723_EEG_NWU_IN', ...
                '250725_EEG_NWU_BN', ...
                '250815_EEG_NWU_PP', ...
                '250819_EEG_NWU_ZL', ...
                '250723_EEG_NWU_BK', ...
                '250912_EEG_NWU_JN', ...
                '250904_OBE_NWU_TI',...
                '250818_Dupi_NMH_JH_2',...
                '250811_Dupi_NMH_TPB_1',...
                '250811_Dupi_NMH_TB_2',...
                '250929_Dupi_NMH_GH_1',...
                '251009_OBE_NWU_CP_1',...
                '251002_Dupi_NMH_AB_1',...
                '251027_Dupi_NMH_DL_1', ...
                '250929_Dupi_NMH_GH_2',...
                '251002_Dupi_NMH_AB_2'...
                '251013_Dupi_NMH_JN_2',...
                '251030_Dupi_NMH_DB_1',...
                '251110_Dupi_NMH_PC_1',...
                '250623_Dupi_NMH_KS_3',...
                '251030_Dupi_NMH_DB_2',...
                '251120_Dupi_NMH_JL_1',...
                '250818_Dupi_NMH_JH_3'};

%participants with new TTL style for more standardized read in:
newList = {'250811_Dupi_NMH_TB_2',...
           '250929_Dupi_NMH_GH_1',...
           '251009_OBE_NWU_CP_1',...
           '251002_Dupi_NMH_AB_1',...
           '251027_Dupi_NMH_DL_1', ...
            '250929_Dupi_NMH_GH_2',...
            '251002_Dupi_NMH_AB_2'...
            '251013_Dupi_NMH_JN_2',...
            '251030_Dupi_NMH_DB_1',...
            '251110_Dupi_NMH_PC_1',...
            '250623_Dupi_NMH_KS_3',...
            '251030_Dupi_NMH_DB_2',...
            '251120_Dupi_NMH_JL_1',...
            '250818_Dupi_NMH_JH_3'};

%there are multiple respiration channels in many recordings
%which one is right for each session: 
rspIDX = [3,3,3,3,1,1,1,1,1,1,3,3,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; 
rspFlip = [-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,1,1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1]; %hard code flip

% addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
% addpath(genpath([codePre 'myFrequentUse']))
% addpath([codePre 'myFrequentUse/export_fig_repo'])

addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))
% addpath([codePre 'fieldtrip-20230118'])
% addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])
addpath([codePre 'ZelanoLabScripts'])

% set(0, 'defaultfigurewindowstyle', 'docked')
% ft_defaults

for sessi = 1:length(sessionIDs)
try
%% custom import for different participants: 
disp(sessi)
%check for pre existing processing: 
if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                sessionIDs{sessi} '_breathingPreProc.mat'], 'file')




    if strcmp(sessionIDs{sessi}, '250811_Dupi_NMH_TPB_1')
        %special handling of JH session 1 because it was recorded in two
        %files
        dat = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                            '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat = dat.curDat; 
        dat2 = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                   '\raw\raw_breathingTasks2/raw_breathingTasks2.mat']);
        dat2 = dat2.curDat; 
        set(0, 'defaultfigurewindowstyle', 'docked')
        
        
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    '250811_Dupi_NMH_TPB_1.csv'];
        behDat = readtable([codePre behDat]);
        
        outDat = struct; 
        outDat.tim = .0005:.0005:300;
        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
    
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat.rawData.trial{1}(idx, :); 
        

       % Find missing samples in the photodiode channel
        nanidx = isnan(photoDiode);
        nNan   = sum(nanidx);
        
        if nNan > 4000
            error('too many missing values!');
        end
        
        % Get raw data for this trial (rows = channels, cols = time)
        rawData = dat.rawData.trial{1};
        
        if nNan > 0
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



        photoDiode = abs(photoDiode); 
        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
        % figure; plot(photoDiode)
        
        TTLs = [112116 709263];
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end

        % % second dataset
        % idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        % photoDiode = dat2.rawData.trial{1}(idx, :); 
        % 
        % photoDiode = abs(photoDiode); 
        % photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
        % figure; plot(photoDiode)
        % 
        % TTLs = [113048, 711003];
        % 
        % TTLs = sort(TTLs); 
        % 
        % for ii = 1:2:length(TTLs)
        %     outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
        %                                       TTLs(ii):TTLs(ii)+599999); 
        %     di = di+1; 
        % end

        %The second dataset is just another round of audio so skip for now

    
    elseif strcmp(sessionIDs{sessi}, '250818_Dupi_NMH_JH_1')
        %special handling of JH session 1 because it was recorded in two
        %files
        dat = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                            '\raw\raw_audioBook/raw_audioBook.mat']);
        dat = dat.curDat; 
        dat2 = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                   '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat2 = dat2.curDat; 
        set(0, 'defaultfigurewindowstyle', 'docked')
        
        
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    '250818_Dupi_NMH_JH_1.csv'];
        behDat = readtable([codePre behDat]);
        
        outDat = struct; 
        
        outDat.data =dat.rawData.trial{1}(:,20000:600000+19999);
        outDat.tim = .0005:.0005:300;
        
        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
        
        
        photoDiodeDat = dat2.rawData.trial{1}; 
        photoDiodeDat = photoDiodeDat(end,:); 
        tim = dat2.rawData.time{1}; 
        %TTLs in sample indices 
        TTLs = find(photoDiodeDat(1:length(photoDiodeDat)-1)<3000 &...
             photoDiodeDat(2:length(photoDiodeDat))>3000);
        
        % figure; plot(photoDiodeDat)
        % xline(TTLs([1, find(diff(TTLs)> 15000), ...
        %                     find(diff(TTLs)> 15000)+1, end]))
        
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                        find(diff(TTLs)> 15000)+1, end]);
        TTLs = sort(TTLs); 
        di = 1; 
        for ii = 1:2:9
            outDat.data(:,:,di+1) = dat2.rawData.trial{1}(:,...
                                               TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end
    
    
    elseif strcmp(sessionIDs{sessi}, '250623_DUPI_NMH_KS_2')
        %SPECIALIZED PROCESSING FOR KS 2
        dat = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat = dat.curDat; 
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '.csv'];
        behDat = readtable([codePre behDat]);
    
     
    
        %the photo diode wasn't working properly, so hard code based on
        %notes
    
        secondBreaks = [210, 510, 620, 920, 1026, 1326,...
                        1380, 1680, 1743, 2043, 2090, 2390];
    
        tim = dat.rawData.time{1}; 
        TTLs = arrayfun(@(x) find(tim>x, 1), secondBreaks); 
    
        outDat = struct; 
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(secondBreaks)/2);
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                             TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end
        outDat.tim = .0005:.0005:300;
        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
    
    elseif strcmp(sessionIDs{sessi}, '250723_EEG_NWU_IN')
        %handle double recording of audio/focus by throwing away the extra
        dat = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat = dat.curDat; 
       
        dat.rawData.trial{1}(:,1:2802720) = []; 
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '.csv'];
        behDat = readtable([codePre behDat]);
        
    
        outDat = struct; 
        outDat.tim = .0005:.0005:300;
        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
    
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat.rawData.trial{1}(idx, :); 
        
    
        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
        % figure; plot(photoDiode)
        TTLs = find(photoDiode(1:length(photoDiode)-1)<1 &...
                                photoDiode(2:length(photoDiode))>1);
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                                find(diff(TTLs)> 15000)+1, end]);
       
        TTLs = sort(TTLs); 
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(TTLs)/2);
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end

   
    elseif strcmp(sessionIDs{sessi}, '250912_EEG_NWU_JN')
        dat = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat = dat.curDat; 
       
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '.csv'];
        behDat = readtable([codePre behDat]);
        
    
        outDat = struct; 
        outDat.tim = .0005:.0005:300;
        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
    
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat.rawData.trial{1}(idx, :); 
        
    
        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
        % figure; plot(photoDiode)
        TTLs = find(photoDiode(1:length(photoDiode)-1)<1 &...
                                photoDiode(2:length(photoDiode))>1);
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                                find(diff(TTLs)> 15000)+1, end]);
       
        TTLs([1,2, 10, 18]) = []; %extra detections for this participant! 

        TTLs = sort(TTLs); 
        % xline(TTLs)
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(TTLs)/2);
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end

    elseif strcmp(sessionIDs{sessi}, '250819_EEG_NWU_ZL')
        dat = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat = dat.curDat; 
       
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '.csv'];
        behDat = readtable([codePre behDat]);
        
    
        outDat = struct; 
        outDat.tim = .0005:.0005:300;
        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
    
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat.rawData.trial{1}(idx, :); 
        
    
        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
        % figure; plot(photoDiode)
        TTLs = find(photoDiode(1:length(photoDiode)-1)<1.1 &...
                                photoDiode(2:length(photoDiode))>1.1);
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                                find(diff(TTLs)> 15000)+1, end]);
       
        

        TTLs = sort(TTLs); 
        % xline(TTLs)
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(TTLs)/2);
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end

     
    elseif sum(cellfun(@(x) strcmp(x, sessionIDs{sessi}), newList))==1 %new standard
        dat = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat = dat.curDat; 
       
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '.csv'];
        behDat = readtable([codePre behDat]);
        
        if strcmp(sessionIDs{sessi} , '251030_Dupi_NMH_DB_2')
            dat.rawData.trial{1}(:,1:1023000) = []; %eliminate initial recording before computer glitch
        end
    
        outDat = struct; 
        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
        outDat.tim = .0005:.0005:300;
    
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat.rawData.trial{1}(idx, :);
        % Find missing samples in the photodiode channel
        nanidx = isnan(photoDiode);
        nNan   = sum(nanidx);
        
        if nNan > 4000
            error('too many missing values!');
        end
        
        % Get raw data for this trial (rows = channels, cols = time)
        rawData = dat.rawData.trial{1};
        
        if nNan > 0
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

        dat.rawData.trial{1} = rawData;
        
        %audio, focused, slow shadow, fast shadow, focus shadow
        cndSeps = [1.1, 1.3, 1.6, 1.7, 1.8, 1.9];

        thresh = prctile(photoDiode, 2);
        TTLs = find(photoDiode(1:length(photoDiode)-1)>thresh &...
                                photoDiode(2:length(photoDiode))<thresh);
        minDist = min(cndSeps) *outDat.fs - 700; 
        TTLs = TTLs([1,  ...
                                find(diff(TTLs)> minDist), end]);
        % TTLs: vector of time stamps (samples or seconds)
        TTLs = TTLs(:);                 % make sure it's a column
        
        % 1) Intervals between consecutive TTLs
        dTTL = diff(TTLs);
        
        % 2) Find the "long" gaps (between blocks)
        %    Here I use a robust outlier rule; tweak 'ThresholdFactor'
        isLongGap = dTTL>outDat.fs*5;
        
        % indices in dTTL of long gaps
        gapIntIdx = find(isLongGap);
        
        % 3) Indices in TTLs:
        idx_before_gap = gapIntIdx;          % last TTL of each block
        idx_after_gap  = gapIntIdx + 1;      % first TTL of next block

        startTTLs = [TTLs(1); TTLs(idx_after_gap)];
        endTTLs = [TTLs(idx_before_gap); TTLs(end)]; 

        blockLens = (endTTLs - startTTLs) ./ outDat.fs;

        TTLs = startTTLs(blockLens>290 & blockLens<310);
        TTLs = sort(TTLs); 

        figure
        plot(photoDiode)
        xline(TTLs)
        title(sessi)
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(TTLs));
        di = 1; 
        for ii = 1:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end


        
        
    else %old STANDARD PROCESSING: 
        dat = load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat = dat.curDat; 
       
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '.csv'];
        behDat = readtable([codePre behDat]);
        
    
        outDat = struct; 
        outDat.tim = .0005:.0005:300;
        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
    
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat.rawData.trial{1}(idx, :); 
        
        photoDiode = abs(photoDiode); 
        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
        % figure; plot(photoDiode)
        
        TTLs = find(photoDiode(1:length(photoDiode)-1)<4 &...
                                photoDiode(2:length(photoDiode))>4);
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                                find(diff(TTLs)> 15000)+1, end]);
       
        TTLs = sort(TTLs); 
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(TTLs)/2);
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end
    
    end
    
    % if strcmp(sessionIDs{sessi},'250811_Dupi_NMH_TPB_1')
    %     outDat.data(:,:,3) = []; 
    % end

    %put all the data into Chan X time with conditions concatenated
    %together
    
    data = zeros(size(outDat.data,1), prod(size(outDat.data,[2,3])));
    for ii = 1:size(outDat.data,3)
        data(:,1+(ii-1)*600000:(ii*600000)) = outDat.data(:,:,ii); 
    end
    
    outDat.data = data; 
    outDat.task = "breathing"; 
    outDat.sessID = sessionIDs{sessi};
    outDat.OGdataDir = [datPre{datPrei(sessi)} sessionIDs{sessi}];
    tmp = dir([datPre{datPrei(sessi)} sessionIDs{sessi}]);
    tmp = tmp(cellfun(@(x) contains(x, '.m'), {tmp.name}));
    tmp = tmp(cellfun(@(x) contains(x, 'LoadData'), {tmp.name}));
    if size(tmp,1) == 1
        outDat.loadFile = tmp.name;
    else 
        error('load file not identified uniquely')
    end
    outDat.preProcScript = 'BreathingTask_makeOutDat.m'; 
    if datPrei(sessi) == 1
        outDat.type = 'Dupi'; 
    elseif datPrei(sessi) == 2
        outDat.type = 'OBE';
    elseif datPrei(sessi) == 3
        outDat.type = 'EEG';
    end



    %make a directory for preprocessed data if it hasn't been made yet
    if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc'], 'dir')
         mkdir([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc']);
    end
    parSave([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_breathingPreProc.mat'], ...
                    outDat)
end


catch 
    disp(['fail for ', num2str(sessi)])
end
end
