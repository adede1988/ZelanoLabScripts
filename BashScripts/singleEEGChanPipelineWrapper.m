
%simulated HPC looping input from shell script

%start = 1;

disp(['attempting file: ' num2str(start)])

%% file path management

%local paths: 

%codePre = 'G:\My Drive\GitHub\';
%datPre = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\';

%HPC paths: 

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

test = cellfun(@(x) length(x)>0, strfind({chanFiles.name}, '_EEG'));
eegFiles = chanFiles(test); 

%cut to breathing task only for now!
test = cellfun(@(x) length(x)>0, strfind({eegFiles.name}, 'breathingTask'));
eegFiles = eegFiles(test); 

%% run the pipeline

% parfor start = 1:length(eegFiles)
%     try
        curChan = eegFiles(start).name; 
        if ~contains(curChan, 'EEG') %check this! 
            subID = split(curChan, '_macro_');
            [~, taskID, ~] = fileparts(subID{2});       
            taskID = regexprep(taskID, '_\d+$', '');   
            subID = subID{1}; 
        else
            subID = split(curChan, '_EEG_'); 
            [~, taskID, ~] = fileparts(subID{2});       
            taskID = regexprep(taskID, '_\d+$', '');   
            subID = subID{1}; 
        end
        
        subFiles = dir([datPre 'CHANDAT']);
        test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, subID)); 
        
        subFiles = subFiles(test);

        test = cellfun(@(x) length(x)>0, strfind({subFiles.name}, taskID)); 
        
        subFiles = subFiles(test);
        
        
        disp(['going for ' subID ' ' eegFiles(start).name] )
        
        if strcmp(taskID, 'breathingTask')
            singleChanEEGPipeline(eegFiles, start, subFiles, datPre); 
        else
            disp('only doing breathing for now, so skip')
        end
    % catch
    %     disp(['failure on ' subID ' ' eegFiles(start).name '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'])
    % end

% end