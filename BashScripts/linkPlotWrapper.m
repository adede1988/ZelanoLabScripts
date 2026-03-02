
%local paths: 

codePre = 'G:\My Drive\GitHub\';
datPre = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\';

%HPC paths: 
% 
% codePre = '/projects/p33197/code/';
% datPre = '/projects/p33197/QuestMirror/';

%% set paths

addpath(genpath([codePre 'ZelanoLabScripts']))
addpath([codePre 'myFrequentUse'])
% addpath([codePre 'fieldtrip-20230118'])
% ft_defaults;

%% initialize chanFiles, dif versions for local v. quest running




% 
datFolder = [datPre 'CHANDAT']; 
chanFiles = dir(datFolder);
test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '.mat'));
chanFiles = chanFiles(test); 

test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '_EEG'));
eegFiles = chanFiles(test); 

%cut to breathing task only for now!
test = cellfun(@(x) length(x)>0, strfind({eegFiles.name}, 'breathingTask'));
eegFiles = eegFiles(test); 
allSubIDs = cell(length(eegFiles),3); 
for ii = 1:length(eegFiles)
    curChan = eegFiles(ii).name;
    subID = split(curChan, '_EEG_'); 
    [~, taskID, ~] = fileparts(subID{2});       
    taskID = regexprep(taskID, '_\d+$', '');   
    subID = subID{1}; 
    allSubIDs{ii,2} = taskID; 
    allSubIDs{ii,1} = subID;
    allSubIDs{ii,3} = [subID '_' taskID];


end



% Create unique integer for each unique (subID, taskID) combination
[comboKey, ~, comboIdx] = unique(allSubIDs(:,3),  'stable');

% comboIdx is your index vector: length == length(eegFiles)
% comboKey is the lookup table: each row is {subID, taskID} for that comboIdx value




%% run the pipeline

for start = 1:max(comboIdx)
    try
        
        subIDX = find(comboIdx == start); 
        subFiles = eegFiles(subIDX);
        
        
        disp(['going for ' allSubIDs{subIDX(1),1} ' ' allSubIDs{subIDX(1),2}] )
        
       %input subFiles, datPre, subID, taskID for macro search
        makeLinkPlots(subFiles, datPre, allSubIDs{subIDX(1),1}, allSubIDs{subIDX(1),2}); 
       
    catch
        disp(['failure on '  allSubIDs{subIDX(1),1} ' ' allSubIDs{subIDX(1),2} '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
    end

end