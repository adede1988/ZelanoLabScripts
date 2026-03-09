
%local paths: 
% 
% codePre = 'G:\My Drive\GitHub\';
% datPre = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\';

%HPC paths: 
% 
codePre = '/projects/p33197/code/';
datPre = '/projects/p33197/QuestMirror/';

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

test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '_macro'));
macFiles = chanFiles(test); 

%cut to breathing task only for now!
test = ~cellfun(@(x) length(x)>0, strfind({macFiles.name}, 'breathingTask'));
macFiles = macFiles(test); 
allSubIDs = cell(length(macFiles),3); 
for ii = 1:length(macFiles)
    curChan = macFiles(ii).name;
    subID = split(curChan, '_macro_'); 
    [~, taskID, ~] = fileparts(subID{2});       
    taskID = regexprep(taskID, '_\d+$', '');   
    subID = subID{1}; 
    allSubIDs{ii,2} = taskID; 
    allSubIDs{ii,1} = subID;
    allSubIDs{ii,3} = [subID '_' taskID];


end



% Create unique integer for each unique (subID, taskID) combination
[comboKey, ~, comboIdx] = unique(allSubIDs(:,3),  'stable');
allSubIDs(:,4) = arrayfun(@(x) x, comboIdx, 'uniformoutput', false); 
% comboIdx is your index vector: length == length(eegFiles)
% comboKey is the lookup table: each row is {subID, taskID} for that comboIdx value




%% run the pipeline

parfor start = 1:max(comboIdx)
    try
        
        subIDX = find(comboIdx == start); 
        subFiles = macFiles(subIDX);
      
        
        disp(['going for ' allSubIDs{subIDX(1),1} ' ' allSubIDs{subIDX(1),2}] )
        
       %input subFiles, datPre, subID, taskID for macro search
        outSum_noEEG_Pipeline(subFiles, datPre, allSubIDs{subIDX(1),1},  allSubIDs{subIDX(1),2})
    catch
        disp(['failure on '  allSubIDs{subIDX(1),1} ' ' allSubIDs{subIDX(1),2} '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
    end

end