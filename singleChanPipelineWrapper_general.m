

%github access token: 
%Feb 20: 
%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


%simulated HPC looping input from shell script

%start = 1;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

% codePre = 'R:\MSS\Johnson_Lab\dtf8829\GitHub\';
% datPre = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\';

%HPC paths: 

codePre = '/projects/p31578/dtf8829/';
datPre = '/projects/p31578/dtf8829/QuestConnect/';

%% set paths

addpath(genpath([codePre 'HpcAccConnectivityProject']))
addpath([codePre 'myFrequentUse'])
addpath([codePre 'fieldtrip-20230118'])
ft_defaults;

%% initialize chanFiles, dif versions for local v. quest running




% 
datFolder = [datPre 'CHANDAT/CHANRAW']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 




%% run the pipeline


curChan = chanFiles(start).name; 
subID = split(curChan, '_'); 
subID = subID{2}; 

subFiles = dir([datPre 'CHANDAT/CHANRAW']);
test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 

subFiles = subFiles(test);



disp(['going for ' subID ' ' chanFiles(start).name] )

singleChanPipeline(chanFiles, start, subFiles, codePre); 



%% put encoding behavioral data onto it
% curSub = 'something'; 
% for ii = 482:length(chanFiles)
%     ii
%     chanDat = load([chanFiles(ii).folder '/' chanFiles(ii).name]).chanDat; 
%     if strcmp(curSub, chanDat.subID)
%         chanDat.encInfo = dat.trialinfo; 
%         parsave([chanFiles(ii).folder '/' chanFiles(ii).name], chanDat)
%     else
%         dat = load(fullfile([chanDat.dataDir '/' chanDat.encDatFn])).data; 
%         chanDat.encInfo = dat.trialinfo; 
%         parsave([chanFiles(ii).folder '/' chanFiles(ii).name], chanDat)
%         curSub = chanDat.subID; 
%     end
%    
% end



% %I need a way of assigning the partner channel 
% partnerChans = struct; 
% pi = 1; 
% tic
% for chan = 1:length(chanFiles)
%    
%     curChan = chanFiles(chan).name; 
%     subID = split(curChan, '_');
% 
%     subID = subID{2}; 
%     
%     subFiles = dir([datPre 'CHANDAT/CHANRAW']);
%     test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 
%     subFiles = subFiles(test);
%     
%     curi = find(cellfun(@(x) strcmp(x, chanFiles(chan).name), {subFiles.name}));
%     if curi <length(subFiles)
%     for chan2 = curi+1:length(subFiles)
%         partnerChans(pi).name = chanFiles(chan).name; 
%         partnerChans(pi).folder = [chanFiles(chan).folder '/CHANRAW']; 
%         partnerChans(pi).name2 = subFiles(chan2).name; 
%         partnerChans(pi).folder2 = subFiles(chan2).folder; 
%         pi = pi+1; 
%         if mod(pi,1000)==0
%             disp(pi)
%             toc
%             tic
%         end
%     end
%     end
% 
% end
% 
% save([codePre 'HpcAccConnectivityProject' '/localChanFiles.mat'], 'partnerChans')
% 
% partnerChans = struct; 
% pi = 1; 
% tic
% for chan = 1:length(chanFiles)
%    
%     curChan = chanFiles(chan).name; 
%     subID = split(curChan, '_');
% 
%     subID = subID{2}; 
%     
%     subFiles = dir([datPre 'CHANDAT/CHANRAW']);
%     test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 
%     subFiles = subFiles(test);
%     
%     curi = find(cellfun(@(x) strcmp(x, chanFiles(chan).name), {subFiles.name}));
%     if curi <length(subFiles)
%     for chan2 = curi+1:length(subFiles)
%         partnerChans(pi).name = chanFiles(chan).name; 
%         partnerChans(pi).folder = [chanFiles(chan).folder '/CHANRAW']; 
%         partnerChans(pi).name2 = subFiles(chan2).name; 
%         partnerChans(pi).folder2 = subFiles(chan2).folder; 
%         pi = pi+1; 
%         if mod(pi,1000)==0
%             disp(pi)
%             toc
%             tic
%         end
%     end
%     end
% 
% end
