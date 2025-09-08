


%custom function to stitch together breathing recordings into one file

clear

codePre = 'G:\My Drive\GitHub\';
datPre = ['R:\Neurology\Zelano_Lab\Lab_Common\' ...
            'Dupi\'];

sessionIDs = {'250818_Dupi_NMH_JH_1', 
               '250623_DUPI_NMH_KS_2',
               '250623_Dupi_NMH_KS_1'};


addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
addpath(genpath([codePre 'myFrequentUse']))
addpath([codePre 'myFrequentUse/export_fig_repo'])

addpath([codePre 'fieldtrip-20230118'])
addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])

ft_defaults

for sessi = 1:length(sessionIDs)

    %% custom import for different participants: 
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


    elseif strcmp(sessionIDs{sessi}, '250623_DUPI_NMH_KS_2')
        %SPECIALIZED PROCESSING FOR KS 2
        dat = load([datPre sessionIDs{sessi} ...
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


    else %STANDARD PROCESSING: 

        dat = load([datPre sessionIDs{sessi} ...
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
        figure; plot(photoDiode)

        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);

        TTLs = find(photoDiode(1:length(photoDiode)-1)<1 &...
                                photoDiode(2:length(photoDiode))>1);
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                                find(diff(TTLs)> 15000)+1, end]);
       
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(TTLs)/2);
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                                    TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end

    end
    
    %% combined flow from here: 

    %at this point there should be an outDat struct with the following: 
    %data: channels X time X conditions
    %tim: 1Xtime vector in seconds
    %behDat: tidy table of behavioral responses to emotion questions
    %labels: 1Xchannels cell array of human readable channel labels
    %CSClist: 1Xchannels cell array of original CSC labels from neuralynx
    %fs: scalar sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull out respiration data: 

    %get out the respiration data
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :, :); 

    %choose which one looks right
    rspDatz = (rspDat - mean(rspDat, [2,3])) ./ std(rspDat, [], [2,3]); 
    thresh = prctile(rspDat, 90, [2,3]); 
    test = arrayfun(@(x) length(find( rspDat(x,1:end-1,3) < thresh(x) & ...
        rspDat(x,2:end,3)>thresh(x))), 1:sum(idx)); 
    [~, idx] = min(test); 
    rspDat = squeeze(rspDat(idx, :, :)); 
    rspDatz = squeeze(rspDatz(idx, :,:)); 

    %check for flipped signal: 
    test = rspDatz(:); 
    if median(test) < -.03
        rspDat = -rspDat; 
    end

    %% get breath summary variables 
    %should be sanity checked by looking at breath traces with onset marks
    %col 1: onset Y value
    %col 2: onset tim
    %col 3: peak Y value
    %col 4: peak tim
    %col 5: end Y value
    %col 6: end tim
    %col 7: length (end tim - onset tim)
    %col 8: amp (peak Y - avg of two ends)
    %col 9: idx of peak in rspSig2
    %col10: exhale peak Y value
    %col11: exhale peak tim
    %col12: condition
    %col13: empty
    %col14: index

    bmObj = breathTemplates4(rspDat(:,1), outDat.fs);
    bmObj(:, 12) = 1; 
    for cndi = 2:size(rspDat, 2)
        bmCur = breathTemplates4(rspDat(:,cndi), outDat.fs);
        bmCur(:,12) = cndi; 
        bmObj = [bmObj; bmCur]; 

    end
   

    for cndi = 1:size(rspDat, 2)
        test = bmObj(bmObj(:,12)==cndi, 2);
        figure
        plot(outDat.tim, rspDat(:, cndi))
        xline(test)

    end

   

    %% get heartbeats: 
        %requires custom heartbeat detection

   
    idx = cellfun(@(x) contains(x, 'ECG'), outDat.labels);
    ECG = outDat.data(idx, :, :); 
    test = reshape(ECG, 3, []); 

    d = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', 5, 'HalfPowerFrequency2', 40, ...
    'SampleRate', outDat.fs);
    
    test = filtfilt(d, test')'; 
    ECG = reshape(test, 3, 600000, []); 

  


    ECGz = (ECG - mean(ECG, [2,3])) ./ std(ECG, [], [2,3]); 
    figure; 
    plot(outDat.tim, ECGz(1,:,1), 'color', 'k')
    hold on 
    plot(outDat.tim,ECGz(2,:,1), 'color', 'red')
    plot(outDat.tim,ECGz(3,:,1), 'color', 'green')

    %working time point: 88.6875
        %% start here!!!!
    test = find(arrayfun(@(x,y,z) x>0 & y>1 & z<-1, ECGz(1,:,1), ...
                                  ECGz(2,:,1), ECGz(3,:,1));
    test = arrayfun(@(x,y,z) x>0 & y>1 & z<-1, ECGz(1,:,1), ...
                                  ECGz(2,:,1), ECGz(3,:,1));

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





