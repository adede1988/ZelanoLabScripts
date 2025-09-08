


%custom function to stitch together breathing recordings into one file

clear

codePre = 'G:\My Drive\GitHub\';
datPre = ['R:\Neurology\Zelano_Lab\Lab_Common\' ...
            'Dupi\'];

sessionIDs = {'250818_Dupi_NMH_JH_1', 
               '250623_DUPI_NMH_KS_2'};


addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
addpath(genpath([codePre 'myFrequentUse']))
addpath([codePre 'myFrequentUse/export_fig_repo'])

addpath([codePre 'fieldtrip-20230118'])
addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])

ft_defaults

for sessi = 1:length(sessionIDs)

    if strcmp(sessionIDs{sessi}, '250818_Dupi_NMH_JH_1')
        %special handling of JH session 1 because it was recorded in two
        %files
        dat = load([datPre sessionIDs{sessi} ...
                            '\raw\raw_audioBook/raw_audioBook.mat']);
        dat = dat.curDat; 
        dat2 = load([datPre sessionIDs{sessi} ...
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
        
        figure; plot(photoDiodeDat)
        xline(TTLs([1, find(diff(TTLs)> 15000), ...
                            find(diff(TTLs)> 15000)+1, end]))
        
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                        find(diff(TTLs)> 15000)+1, end]);
        di = 1; 
        for ii = 1:2:9
            outDat.data(:,:,di+1) = dat2.rawData.trial{1}(:,...
                                                    TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end


    else
        dat = load([datPre sessionIDs{sessi} ...
                       '\raw\raw_breathingTasks/raw_breathingTasks.mat']);
        dat = dat.curDat; 
        behDat = ['closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '.csv'];
        behDat = readtable([codePre behDat]);

        outDat.behDat = behDat; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 

        photoDiodeDat = dat.rawData.trial{1}; 
        photoDiodeDat = photoDiodeDat(end,:); 
        tim = dat.rawData.time{1}; 
        %TTLs in sample indices 
        TTLs = find(photoDiodeDat(1:length(photoDiodeDat)-1)<3000 &...
             photoDiodeDat(2:length(photoDiodeDat))>3000);
        
        figure; plot(photoDiodeDat)


    end





end


%Figure generation for Tate: 

% tLen = dat.rawData.time{1}; 
% winSize = 20000;
% breaks = 1:winSize:length(tLen);
% 
% for ii = 1:length(breaks)-1
%     figure('visible', false, 'position', [0,0,1600,1000])
%     hold on
% 
%     for chan = 1:6
%         plot(tLen(breaks(ii):breaks(ii+1)), ...
%             dat.rawData.trial{1}(32+chan, breaks(ii):breaks(ii+1)) + 200*chan, ...
%             'color', 'blue', 'linewidth', 1.5)
%     end
% 
%     for chan = 1:5:32
%         plot(tLen(breaks(ii):breaks(ii+1)), ...
%             dat.rawData.trial{1}(chan, breaks(ii):breaks(ii+1)) + ...
%             1400 + 200*(chan-1)/5, ...
%             'color', 'green', 'linewidth', 1.5)
%     end
% 
%     yticks([200:200:2600])
%     yticklabels({'macro1', 'macro2', 'macro3', 'macro4', 'macro5', 'macro6', ...
%         dat.outLabs{1:5:32}})
% 
%     title(['file: ' num2str(ii)])
%     exampDat = struct;
%     exampDat.data = dat.rawData.trial{1}([1:5:32, 33:38], breaks(ii):breaks(ii+1));
%     exampDat.tim = tLen(breaks(ii):breaks(ii+1));
%     exampDat.labels = {'macro1', 'macro2', 'macro3', 'macro4', 'macro5', 'macro6', ...
%         dat.outLabs{1:5:32}};
% 
%     save(['G:\My Drive\GitHub\ZelanoLabScripts\noiseFigs\' ...
%         'NoiseExample_' num2str(ii) '.mat'], 'exampDat')
%     export_fig(['G:\My Drive\GitHub\ZelanoLabScripts\noiseFigs\' ...
%         'NoiseExample_' num2str(ii) '.jpg'], '-r300')
% end





