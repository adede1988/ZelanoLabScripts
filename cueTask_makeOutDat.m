clear

behDatPath = 'R:\Neurology\Zelano_Lab\Lab_Common\OBE_task_backup\tasks\olf_cuetask\results\olf_cuetask_results';
codePre = 'C:\Users\Adam\Documents\GitHub\';

% addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
% addpath(genpath([codePre 'myFrequentUse']))
addpath(genpath([codePre 'ZelanoLabScripts']))
% addpath([codePre 'myFrequentUse/export_fig_repo'])
addpath(genpath('C:\Users\Adam\Documents\eeglab2026.0.0'))

% addpath([codePre 'fieldtrip-20230118'])
% addpath([codePre 'emotionDecoding'])
% addpath([codePre 'slowBreathing'])

set(0, 'defaultfigurewindowstyle', 'docked')
% ft_defaults

%% 

cfg        = applyParams('cueTask','makeOutDat');
sessionIDs = cfg.sessionIDs;
datPre     = cfg.datPre;
datPrei    = cfg.datPrei;
newSet     = cfg.newIDs;
rspIDX     = cfg.rspIDX;
rspFlip    = cfg.rspFlip;

for sessi = 1:length(sessionIDs)
    disp(['Working on: ' sessionIDs{sessi}])

    

    if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_cueTaskPreProc.mat'], 'file')
%% new set 
if sum(cellfun(@(x) strcmp(sessionIDs{sessi}, x), newSet))==1

    datFolders = dir([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\']);

    idx = cellfun(@(x) contains(x, 'raw_cueTaskOdor'), {datFolders.name});
    idx = find(idx); 
    
    %new behavioral data location 
    dat1 = load([datFolders(idx(1)).folder '\' ...
                     datFolders(idx(1)).name  '\' ...
                     datFolders(idx(1)).name  '.mat']);
    dat1 = dat1.curDat; 

       

          %find the behavioral files: 
        behDir = dir([behDatPath, '\' sessionIDs{sessi}]);
       
       
        idx = cellfun(@(x) ~contains(x, 'mat'), {behDir.name});
        behDir(idx) = [];
        
        if length(behDir) ~= 2
            error('wrong number of behavior files available!')
        end
        %load the behavioral files
        behDat1 = [behDir(1).folder '\' behDir(1).name];
        behDat1 = load(behDat1);
        behDat2 = [behDir(2).folder '\' behDir(2).name];
        behDat2 = load(behDat2);
        
        behDat1 = outMat_to_table(behDat1.outMat, behDat1.datalabel);
        behDat2 = outMat_to_table(behDat2.outMat, behDat2.datalabel);

        %check if there is a value greater than 10 in the odor column
        %recode if necessary
        if max(behDat1.odor) == 11 || min(behDat1.odor)==2
            behDat1.odor = behDat1.odor - 1; 
            behDat2.odor = behDat2.odor - 1; 
        end

         %start making the oudDat struct
        outDat = struct; 
        outDat.labels = dat1.outLabs;
        outDat.CSClist = dat1.ncslabels; 
        outDat.fs = dat1.rawData.fsample; 
        outDat.sessID = sessionIDs{sessi}; 

        
        %process the TTL pulses for dat 1: 
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat1.rawData.trial{1}(idx, :); 

         % Find missing samples in the photodiode channel
        nanidx = isnan(photoDiode);
        nNan   = sum(nanidx);
        
        if nNan > 4000
            error('too many missing values!');
        end
        
        % Get raw data for this trial (rows = channels, cols = time)
        rawData = dat1.rawData.trial{1};
        
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

        dat1.rawData.trial{1} = rawData;

        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
    

        downs = find(photoDiode(1:end-1) > -1.5 & photoDiode(2:end)<-1.5);
        ups = find(photoDiode(1:end-1) < -1.5 & photoDiode(2:end)>-1.5);
        difVals = ups - downs; %difVals is length of TTL pulses
        downs(difVals>1500) = []; 
        ups(difVals>1500) = []; 
        difVals(difVals>1500) = []; 

        downs(difVals<300) = []; 
        ups(difVals<300) = []; 
        difVals(difVals<300) = []; 

        starti = []; 
        for di = 5:length(downs)
            if downs(di) - downs(di-4) < 3500
                starti = [starti di]; 
            end
        end

        %find two key starts
        if length(starti)>1
            spliti = find(diff(starti)>10);
            starti1 = starti(spliti); 
            starti(1:spliti) = []; 
            if length(starti)>1
                spliti = find(diff(starti)>10);
                if isempty(spliti)
                    starti = [starti1 starti(end)]; 
                else
                    error('unknown TTL structure')
                end
            else 
                error('only one run?')
            end
        end

        %eliminate pre start 1 TTLs
        downs(1:starti(1)) = []; 
        difVals(1:starti(1)) = []; 

        %find last TTL of run 1
        tmp = diff(downs);
        [~, splitIDX] = max(tmp);

        %eliminate pre start 2 TTLs
        downs(splitIDX-1:starti(2)-starti(1)) = []; 
        difVals(splitIDX-1:starti(2)-starti(1)) = []; 

        %eliminate trailing TTL
        downs(end-1:end) = []; 
        difVals(end-1:end) = []; 

        %count missed responses to know how many TTLs there should be
        
        s1 = sum(arrayfun(@(x) x==999, behDat1.response));
        s2 = sum(arrayfun(@(x) x==999, behDat2.response));
       
        if length(downs) ~= (120 - s1 - s2)
            error('wrong number of TTLs')
        end
        
        cueTTLs = downs(difVals < 850);
        sniffTTLs = downs(difVals > 850 & difVals < 1250);
        responseTTLs = downs(difVals > 1250); 
    
        figure
        plot(photoDiode)
        xline(cueTTLs, 'color', 'magenta', 'linewidth', 2)
        xline(sniffTTLs, 'color', 'green', 'linewidth', 2)
        xline(responseTTLs, 'color', 'k', 'linewidth',2)


        behDat = [behDat1; behDat2];
        %store all TTLs into one matrix: 
        %col 1: trialStarts
        %col 2: sniff 
        %col 3: buttonPress
       
        TTLs = nan(40, 3); 
        ri = 1; 
        for tt = 1:40
            if isempty(behDat.response_str{tt})
                TTLs(tt,1) = cueTTLs(tt); 
                TTLs(tt,2) = sniffTTLs(tt); 
            else
                TTLs(tt,1) = cueTTLs(tt); 
                TTLs(tt,2) = sniffTTLs(tt); 
                TTLs(tt,3) = responseTTLs(ri);
                ri = ri+1; 
            end
        end
        TTLs =round(TTLs ./ 4);

        outDat.behDat = behDat; 
        outDat.data = dat1.rawData.trial{1};





else
   %% not part of the new set
    %find behavioral folder: 
    subFolders = dir([datPre{datPrei(sessi)} sessionIDs{sessi}]);
    subFolders = subFolders([subFolders.isdir]);
    idx = cellfun(@(x) contains(x, 'ehavior'), {subFolders.name});
    idx = find(idx); 
    if length(idx) ~= 1
        error('behavioral data folder not uniquely identified!')
    end
    behFold = subFolders(idx);
    datFolders = dir([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                       '\raw\']);

    idx = cellfun(@(x) contains(x, 'raw_cueTaskOdor'), {datFolders.name});
    idx = find(idx); 
    
%% double file participants: 
    if length(idx) == 2
        %process when two runs need to be stitched together: 

        %load in both files: 
        dat1 = load([datFolders(idx(1)).folder '\' ...
                     datFolders(idx(1)).name  '\' ...
                     datFolders(idx(1)).name  '.mat']);
        dat1 = dat1.curDat; 
        dat2 = load([datFolders(idx(2)).folder '\' ...
                     datFolders(idx(2)).name  '\' ...
                     datFolders(idx(2)).name  '.mat']);
        dat2 = dat2.curDat; 

        %behavioral data loading
        switch sessionIDs{sessi}
            case '230611_OBE_NMH_AZ'
                %find the behavioral files: 
                behDir = dir([behFold.folder, '\' ...
                              behFold.name '\olf_cue']);
                %eliminate echem and imagine files
                idx = cellfun(@(x) contains(x, 'echem'), {behDir.name});
                idx = find(idx);
                behDir(idx) = []; 
                idx = cellfun(@(x) contains(x, 'imagine'), {behDir.name});
                idx = find(idx);
                behDir(idx) = [];
                idx = cellfun(@(x) contains(x, 'txt'), {behDir.name});
                idx = find(idx);
                behDir(idx) = [];
                idx = cellfun(@(x) contains(x, 'mat'), {behDir.name});
                idx = find(idx);
                if length(idx) ~= 3
                    error('wrong number of behavior files available!')
                end
                %load the behavioral files
                behDat1 = [behDir(idx(1)).folder '\' behDir(idx(1)).name];
                behDat1 = load(behDat1);
                behDat2 = [behDir(idx(2)).folder '\' behDir(idx(2)).name];
                behDat2 = load(behDat2);
                labs = [behDir(idx(3)).folder '\' behDir(idx(3)).name];
                labs = load(labs);

                behDat1 = outMat_to_table(behDat1.outMat, labs.datalabel);
                behDat2 = outMat_to_table(behDat2.outMat, labs.datalabel);
    
            case '250310_OBE_NMH_FS'
                 %find the behavioral files: 
                behDir = dir([behFold.folder, '\' ...
                              behFold.name '\olf_cue']);
                %eliminate echem and imagine files
                idx = cellfun(@(x) contains(x, 'echem'), {behDir.name});
                idx = find(idx);
                behDir(idx) = []; 
                idx = cellfun(@(x) contains(x, 'imagine'), {behDir.name});
                idx = find(idx);
                behDir(idx) = [];
                idx = cellfun(@(x) contains(x, 'txt'), {behDir.name});
                idx = find(idx);
                behDir(idx) = [];
                idx = cellfun(@(x) contains(x, 'mat'), {behDir.name});
                idx = find(idx);
                if length(idx) ~= 2
                    error('wrong number of behavior files available!')
                end
                %load the behavioral files
                behDat1 = [behDir(idx(1)).folder '\' behDir(idx(1)).name];
                behDat1 = load(behDat1);
                behDat2 = [behDir(idx(2)).folder '\' behDir(idx(2)).name];
                behDat2 = load(behDat2);
               

                behDat1 = outMat_to_table(behDat1.outMat, ...
                                          behDat1.datalabel);
                behDat2 = outMat_to_table(behDat2.outMat, ...
                                          behDat2.datalabel);

            otherwise

                %find the behavioral files: 
                behDir = dir([behFold.folder, '\' behFold.name '\olf_cue']);
                %eliminate echem and imagine files
                idx = cellfun(@(x) contains(x, 'echem'), {behDir.name});
                idx = find(idx);
                behDir(idx) = []; 
                idx = cellfun(@(x) contains(x, 'imagine'), {behDir.name});
                idx = find(idx);
                behDir(idx) = [];
                idx = cellfun(@(x) contains(x, 'visual'), {behDir.name});
                idx = find(idx);
                behDir(idx) = [];
                idx = cellfun(@(x) contains(x, 'mat'), {behDir.name});
                idx = find(idx);
                behDir(idx) = [];
                idx = cellfun(@(x) contains(x, 'txt'), {behDir.name});
                idx = find(idx);
                if length(idx) ~= 2
                    error('wrong number of behavior files available!')
                end
                %load the behavioral files
                behDat1 = [behDir(idx(1)).folder '\' behDir(idx(1)).name];
                behDat1 = readtable(behDat1);
                behDat2 = [behDir(idx(2)).folder '\' behDir(idx(2)).name];
                behDat2 = readtable(behDat2);
        end
        %check if there is a value greater than 10 in the odor column
        %recode if necessary
        if max(behDat1.odor) == 11 || min(behDat1.odor)==2
            behDat1.odor = behDat1.odor - 1; 
            behDat2.odor = behDat2.odor - 1; 
        end



        %start making the oudDat struct
        outDat = struct; 
        outDat.labels = dat1.outLabs;
        outDat.CSClist = dat1.ncslabels; 
        outDat.fs = dat1.rawData.fsample; 
        outDat.sessID = sessionIDs{sessi}; 

        
        %process the TTL pulses for dat 1: 
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat1.rawData.trial{1}(idx, :); 

        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
    

        downs = find(photoDiode(1:end-1) > -1.5 & photoDiode(2:end)<-1.5);
        ups = find(photoDiode(1:end-1) < -1.5 & photoDiode(2:end)>-1.5);
        difVals = ups - downs; %difVals is length of TTL pulses
        downs(difVals>1500) = []; 
        ups(difVals>1500) = []; 
        difVals(difVals>1500) = []; 

        downs(difVals<300) = []; 
        ups(difVals<300) = []; 
        difVals(difVals<300) = []; 

        switch sessionIDs{sessi}
            case '230611_OBE_NMH_AZ'
       
            

            case  '241017_OBE_NMH_AS'
                starti = 1; 
                for di = 5:length(downs)
                    if downs(di) - downs(di-4) < 3500
                        starti = di; 
                    end
                end
                downs(1:starti) = []; 
                difVals(1:starti) = []; 
                downs(end-1:end) = []; 
                difVals(end-1:end) = []; 
        
                if length(downs) ~= 40
                    error('wrong number of TTLs')
                end
            otherwise 
                starti = 1; 
                for di = 5:length(downs)
                    if downs(di) - downs(di-4) < 3500
                        starti = di; 
                    end
                end
                downs(1:starti) = []; 
                difVals(1:starti) = []; 
                downs(end-1:end) = []; 
                difVals(end-1:end) = []; 
        
                if length(downs) ~= 60 - sum(arrayfun(@(x) x==999, ...
                        behDat1.response))
                    error('wrong number of TTLs')
                end
        end

        cueTTLs = downs(difVals < 850);
        sniffTTLs = downs(difVals > 850 & difVals < 1250);
        responseTTLs = downs(difVals > 1250); 
    
        figure
        plot(photoDiode)
        xline(cueTTLs, 'color', 'magenta', 'linewidth', 2)
        xline(sniffTTLs, 'color', 'green', 'linewidth', 2)
        
        if isempty(responseTTLs)
            TTLs = nan(20, 3); 
            TTLs(:,1) = cueTTLs; 
            TTLs(:,2) = sniffTTLs; 
                TTLs =round(TTLs ./ 4);
        else
            try
                xline(responseTTLs, 'color', 'k', 'linewidth',2)
                %store all TTLs into one matrix: 
                %col 1: trialStarts
                %col 2: sniff 
                %col 3: buttonPress
                
                TTLs = nan(20, 3); 
                TTLs(:,1) = cueTTLs; 
                TTLs(:,2) = sniffTTLs; 
                TTLs(:,3) = responseTTLs;  
                TTLs =round(TTLs ./ 4);
            catch 
                TTLs = nan(20, 3); 
                ri = 1; 
                for tt = 1:size(behDat1,1)
                    if isempty(behDat1.response_str{tt})
                        TTLs(tt,1) = cueTTLs(tt); 
                        TTLs(tt,2) = sniffTTLs(tt); 
                    else
                        TTLs(tt,1) = cueTTLs(tt); 
                        TTLs(tt,2) = sniffTTLs(tt); 
                        TTLs(tt,3) = responseTTLs(ri);
                        ri = ri+1; 
                    end
                end
                TTLs =round(TTLs ./ 4);
                warning('missed response TTLs')
                
            end
        end

        %process the TTL pulses for dat 2: 
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat2.rawData.trial{1}(idx, :); 

        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
    

        downs = find(photoDiode(1:end-1) > -1.5 & photoDiode(2:end)<-1.5);
        ups = find(photoDiode(1:end-1) < -1.5 & photoDiode(2:end)>-1.5);
        difVals = ups - downs; %difVals is length of TTL pulses
        downs(difVals>1500) = []; 
        ups(difVals>1500) = []; 
        difVals(difVals>1500) = []; 

        downs(difVals<300) = []; 
        ups(difVals<300) = []; 
        difVals(difVals<300) = []; 

         switch sessionIDs{sessi}
            case '230611_OBE_NMH_AZ'
                downs(1:3) = []; 
                difVals(1:3) = []; 
            

            case  '241017_OBE_NMH_AS'
                starti = 1; 
                for di = 5:length(downs)
                    if downs(di) - downs(di-4) < 3500
                        starti = di; 
                    end
                end
                downs(1:starti) = []; 
                difVals(1:starti) = []; 
                downs(end-1:end) = []; 
                difVals(end-1:end) = []; 
        
                if length(downs) ~= 40
                    error('wrong number of TTLs')
                end
            otherwise 
                starti = 1; 
                for di = 5:length(downs)
                    if downs(di) - downs(di-4) < 3500
                        starti = di; 
                    end
                end
                downs(1:starti) = []; 
                difVals(1:starti) = []; 
                downs(end-1:end) = []; 
                difVals(end-1:end) = []; 
        
                if length(downs) ~= 60 - sum(arrayfun(@(x) x==999,...
                                behDat2.response))
                    error('wrong number of TTLs')
                end
         end

        cueTTLs = downs(difVals < 850);
        sniffTTLs = downs(difVals > 850 & difVals < 1250);
        responseTTLs = downs(difVals > 1250); 
    
         figure
        plot(photoDiode)
        xline(cueTTLs, 'color', 'magenta', 'linewidth', 2)
        xline(sniffTTLs, 'color', 'green', 'linewidth', 2)
          if isempty(responseTTLs)
            TTLs2 = nan(20, 3); 
            TTLs2(:,1) = cueTTLs; 
            TTLs2(:,2) = sniffTTLs; 
                TTLs2 =round(TTLs2 ./ 4);
        else
            try
                xline(responseTTLs, 'color', 'k', 'linewidth',2)
                %store all TTLs into one matrix: 
                %col 1: trialStarts
                %col 2: sniff 
                %col 3: buttonPress
               
                TTLs2 = nan(20, 3); 
                TTLs2(:,1) = cueTTLs; 
                TTLs2(:,2) = sniffTTLs; 
                TTLs2(:,3) = responseTTLs;  
                TTLs2 =round(TTLs2 ./ 4);
            catch 
                TTLs2 = nan(20, 3); 
                ri = 1; 
                for tt = 1:size(behDat1,1)
                    if isempty(behDat1.response_str{tt})
                        TTLs2(tt,1) = cueTTLs(tt); 
                        TTLs2(tt,2) = sniffTTLs(tt); 
                    else
                        TTLs2(tt,1) = cueTTLs(tt); 
                        TTLs2(tt,2) = sniffTTLs(tt); 
                        TTLs2(tt,3) = responseTTLs(ri);
                        ri = ri+1; 
                    end
                end
                TTLs2 =round(TTLs2 ./ 4);
                warning('missed response TTLs')
                
            end
          end
        

        %combine everything! 
        L1 = size(dat1.rawData.trial{1},2);
        TTLs2 = TTLs2 + round(L1/4); 

        TTLs = [TTLs; TTLs2]; 
        comboDat = [dat1.rawData.trial{1}, ...
                    dat2.rawData.trial{1}];
        behDat = [behDat1; behDat2]; 

        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = comboDat(idx, :); 

        figure
        plot(photoDiode)
        xline(TTLs(:,1).*4, 'color', 'magenta')
        xline(TTLs(:,2).*4, 'color', 'green')
        xline(TTLs(:,3).*4, 'color', 'k')
        
        outDat.behDat = behDat; 
        outDat.data = comboDat;
        
%% Single file
    else   
        dat1 = load([datFolders(idx(1)).folder '\' ...
                     datFolders(idx(1)).name  '\' ...
                     datFolders(idx(1)).name  '.mat']);
        dat1 = dat1.curDat; 

       

          %find the behavioral files: 
        behDir = dir([behFold.folder, '\' behFold.name '\OdorCueTask']);
        %eliminate echem and imagine files
        idx = cellfun(@(x) contains(x, 'echem'), {behDir.name});
        idx = find(idx);
        behDir(idx) = []; 
        idx = cellfun(@(x) contains(x, 'imagine'), {behDir.name});
        idx = find(idx);
        behDir(idx) = [];
        idx = cellfun(@(x) contains(x, 'txt'), {behDir.name});
        idx = find(idx);
        behDir(idx) = [];
        idx = cellfun(@(x) contains(x, 'mat'), {behDir.name});
        idx = find(idx);
        if length(idx) ~= 2
            error('wrong number of behavior files available!')
        end
        %load the behavioral files
        behDat1 = [behDir(idx(1)).folder '\' behDir(idx(1)).name];
        behDat1 = load(behDat1);
        behDat2 = [behDir(idx(2)).folder '\' behDir(idx(2)).name];
        behDat2 = load(behDat2);
        
        behDat1 = outMat_to_table(behDat1.outMat, behDat1.datalabel);
        behDat2 = outMat_to_table(behDat2.outMat, behDat2.datalabel);

        %check if there is a value greater than 10 in the odor column
        %recode if necessary
        if max(behDat1.odor) == 11 || min(behDat1.odor)==2
            behDat1.odor = behDat1.odor - 1; 
            behDat2.odor = behDat2.odor - 1; 
        end

         %start making the oudDat struct
        outDat = struct; 
        outDat.labels = dat1.outLabs;
        outDat.CSClist = dat1.ncslabels; 
        outDat.fs = dat1.rawData.fsample; 
        outDat.sessID = sessionIDs{sessi}; 

        
        %process the TTL pulses for dat 1: 
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat1.rawData.trial{1}(idx, :); 


         % Find missing samples in the photodiode channel
        nanidx = isnan(photoDiode);
        nNan   = sum(nanidx);
        
        if nNan > 4000
            error('too many missing values!');
        end
        
        % Get raw data for this trial (rows = channels, cols = time)
        rawData = dat1.rawData.trial{1};
        
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

        dat1.rawData.trial{1} = rawData;






        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
    

        downs = find(photoDiode(1:end-1) > -1.5 & photoDiode(2:end)<-1.5);
        ups = find(photoDiode(1:end-1) < -1.5 & photoDiode(2:end)>-1.5);
        difVals = ups - downs; %difVals is length of TTL pulses
        downs(difVals>1500) = []; 
        ups(difVals>1500) = []; 
        difVals(difVals>1500) = []; 

        downs(difVals<300) = []; 
        ups(difVals<300) = []; 
        difVals(difVals<300) = []; 

        starti = []; 
        for di = 5:length(downs)
            if downs(di) - downs(di-4) < 3500
                starti = [starti di]; 
            end
        end

        %find two key starts
        if length(starti)>1
            spliti = find(diff(starti)>10);
            starti1 = starti(spliti); 
            starti(1:spliti) = []; 
            if length(starti)>1
                spliti = find(diff(starti)>10);
                if isempty(spliti)
                    starti = [starti1 starti(end)]; 
                else
                    error('unknown TTL structure')
                end
            else 
                error('only one run?')
            end
        end

        %eliminate pre start 1 TTLs
        downs(1:starti(1)) = []; 
        difVals(1:starti(1)) = []; 

        %find last TTL of run 1
        tmp = diff(downs);
        [~, splitIDX] = max(tmp);

        %eliminate pre start 2 TTLs
        downs(splitIDX-1:starti(2)-starti(1)) = []; 
        difVals(splitIDX-1:starti(2)-starti(1)) = []; 

        %eliminate trailing TTL
        downs(end-1:end) = []; 
        difVals(end-1:end) = []; 

        %count missed responses to know how many TTLs there should be
        
        s1 = sum(arrayfun(@(x) x==999, behDat1.response));
        s2 = sum(arrayfun(@(x) x==999, behDat2.response));
       
        if length(downs) ~= (120 - s1 - s2)
            error('wrong number of TTLs')
        end
        
        cueTTLs = downs(difVals < 850);
        sniffTTLs = downs(difVals > 850 & difVals < 1250);
        responseTTLs = downs(difVals > 1250); 
    
        figure
        plot(photoDiode)
        xline(cueTTLs, 'color', 'magenta', 'linewidth', 2)
        xline(sniffTTLs, 'color', 'green', 'linewidth', 2)
        xline(responseTTLs, 'color', 'k', 'linewidth',2)


        behDat = [behDat1; behDat2];
        %store all TTLs into one matrix: 
        %col 1: trialStarts
        %col 2: sniff 
        %col 3: buttonPress
       
        TTLs = nan(40, 3); 
        ri = 1; 
        for tt = 1:40
            if isempty(behDat.response_str{tt})
                TTLs(tt,1) = cueTTLs(tt); 
                TTLs(tt,2) = sniffTTLs(tt); 
            else
                TTLs(tt,1) = cueTTLs(tt); 
                TTLs(tt,2) = sniffTTLs(tt); 
                TTLs(tt,3) = responseTTLs(ri);
                ri = ri+1; 
            end
        end
        TTLs =round(TTLs ./ 4);

        outDat.behDat = behDat; 
        outDat.data = dat1.rawData.trial{1};
    end
end
    %process behavioral data into hit/miss/cr/fa
    outDat.behDat.type = repmat("hit", 40,1); 
    idx = find(outDat.behDat.cue==outDat.behDat.odor & ...
                       strcmp(outDat.behDat.response_str,'No'));
    outDat.behDat.type(idx) = repmat("miss", length(idx),1); 
    idx = find(outDat.behDat.cue~=outDat.behDat.odor & ...
                       strcmp(outDat.behDat.response_str,'No'));
    outDat.behDat.type(idx) = repmat("cr", length(idx),1); 
    idx = find(outDat.behDat.cue~=outDat.behDat.odor & ...
                       strcmp(outDat.behDat.response_str,'Yes'));
    outDat.behDat.type(idx) = repmat("fa", length(idx),1);  
    idx = find(cellfun(@(x) isempty(x), outDat.behDat.response_str));
    outDat.behDat.type(idx) = repmat("skip", length(idx),1);  
  
    
    outDat.task = "cueTask"; 
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
    outDat.preProcScript = 'cueTaskPreProc.m'; 
    if datPrei(sessi) == 1
        outDat.type = 'Dupi'; 
    elseif datPrei(sessi) == 2
        outDat.type = 'OBE';
    end
    outDat.TTL = table; 
    outDat.TTL.trialStart = TTLs(:,1);  
    outDat.TTL.response = TTLs(:,3); 
    outDat.TTL.sniff = TTLs(:,2);


    save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_cueTaskPreProc.mat'], ...
                    'outDat', "-v7.3")



    clear behDat behDat1 behDat2 behDir behFold cueTTLs dat1 datalabel ...
        datFolders di difVals downs idx outMat photoDiode responseTTLs ...
        ri s1 s2 sniffTTLs spliti splitIDX starti starti1 subFolders ...
        tmp tt ups TTLs2 L1 comboDat dat2 labs
    else
        disp(['already done on: ' sessionIDs{sessi}])
    end

end

