


%custom function to stitch together breathing recordings into one file

clear

codePre = 'G:\My Drive\GitHub\';
datPre = { 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ... 
           'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\',...
           'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

%prefix index for data folder: 
datPrei = [1,1,1,2,3,3,3,3,3,3,2,1,1]; 

sessionIDs = {'250818_Dupi_NMH_JH_1', ...%preprocessed
               '250623_DUPI_NMH_KS_2',...%preprocessed
               '250623_Dupi_NMH_KS_1',...%preprocessed
               '250908_OBE_NWU_AS', ...%preprocessed
                '250723_EEG_NWU_IN', ...%preprocessed
                '250725_EEG_NWU_BN', ...%preprocessed
                '250815_EEG_NWU_PP', ...%preprocessed
                '250819_EEG_NWU_ZL', ...%preprocessed
                '250723_EEG_NWU_BK', ...%preprocessed
                '250912_EEG_NWU_JN', ...%preprocessed
                '250904_OBE_NWU_TI',...%preprocessed
                '250818_Dupi_NMH_JH_2',...%preprocessed
                '250811_Dupi_NMH_TPB_1'};%preprocessed

%there are multiple respiration channels in many recordings
%which one is right for each session: 
rspIDX = [3,3,3,3,1,1,1,1,1,1,3,3,3]; 
rspFlip = [-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,1,1]; %hard code flip

addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
addpath(genpath([codePre 'myFrequentUse']))
addpath([codePre 'myFrequentUse/export_fig_repo'])

addpath([codePre 'fieldtrip-20230118'])
addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])

set(0, 'defaultfigurewindowstyle', 'docked')
ft_defaults

for sessi = 1:length(sessionIDs)

%% custom import for different participants: 

%check for pre existing processing: 
if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' 'TESTTESTTESTTESTTEST' ...
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
        
        photoDiode = abs(photoDiode); 
        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
        figure; plot(photoDiode)
        
        TTLs = find(photoDiode(1:length(photoDiode)-1)>4 &...
                                photoDiode(2:length(photoDiode))>4);
        keep = zeros(size(TTLs));
        prior = -Inf; 
        for ii = 1:length(TTLs)
            if TTLs(ii) - prior > 100000
                keep(ii) = 1; 
            end
                prior = TTLs(ii); 
        end
        idx = find(keep);
        idx = [idx, idx - 1]; 
        idx(idx == 0) = []; 
        idx = sort(idx); 
        idx(end) = []; 
        xline(TTLs(idx))
        TTLs = TTLs(idx);
       
        TTLs = sort(TTLs); 
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(TTLs)/2);
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end

        % second dataset
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat2.rawData.trial{1}(idx, :); 

        photoDiode = abs(photoDiode); 
        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
        figure; plot(photoDiode)

        TTLs = [113048, 711003];

        TTLs = sort(TTLs); 

        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end

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
        
        figure; plot(photoDiodeDat)
        xline(TTLs([1, find(diff(TTLs)> 15000), ...
                            find(diff(TTLs)> 15000)+1, end]))
        
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
        figure; plot(photoDiode)
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
        figure; plot(photoDiode)
        TTLs = find(photoDiode(1:length(photoDiode)-1)<1 &...
                                photoDiode(2:length(photoDiode))>1);
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                                find(diff(TTLs)> 15000)+1, end]);
       
        TTLs([1,2, 10, 18]) = []; %extra detections for this participant! 

        TTLs = sort(TTLs); 
        xline(TTLs)
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
        figure; plot(photoDiode)
        TTLs = find(photoDiode(1:length(photoDiode)-1)<1.1 &...
                                photoDiode(2:length(photoDiode))>1.1);
        TTLs = TTLs([1, find(diff(TTLs)> 15000), ...
                                find(diff(TTLs)> 15000)+1, end]);
       
        

        TTLs = sort(TTLs); 
        xline(TTLs)
        outDat.data = zeros(size(dat.rawData.trial{1}, 1), 600000, ...
                    length(TTLs)/2);
        di = 1; 
        for ii = 1:2:length(TTLs)
            outDat.data(:,:,di) = dat.rawData.trial{1}(:,...
                                              TTLs(ii):TTLs(ii)+599999); 
            di = di+1; 
        end

     
    else %STANDARD PROCESSING: 
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
        figure; plot(photoDiode)
        
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
    
    %make a directory for preprocessed data if it hasn't been made yet
    if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc'], 'dir')
         mkdir([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc']);
    end
    % save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
    %                 sessionIDs{sessi} '_breathingPreProc.mat'], ...
    %                 'outDat', "-v7.3")
%if there is pre existing processing, then load it: 
else
    outDat =  load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
                          '\preProc\' sessionIDs{sessi} ...
                          '_breathingPreProc.mat']).outDat; 

end
    
    %% combined flow from here: 

    %at this point there should be an outDat struct with the following: 
    %data: channels X time X conditions
    %tim: 1Xtime vector in seconds
    %behDat: tidy table of behavioral responses to emotion questions
    %labels: 1Xchannels cell array of human readable channel labels
    %CSClist: 1Xchannels cell array of original CSC labels from neuralynx
    %fs: scalar sampling rate

outDat.sessID = sessionIDs{sessi}; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% downsample and bandpass: 
outDat = downsample_data(outDat, 500);


%% data cleaning
[C,T,N] = size(outDat.data);
X = reshape(outDat.data, C, T*N);  

try
    idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
    macroDat = X(idx, :); 
    macOut = zeros([5, size(macroDat, 2)]);
    %do bipolar rereferencing 
    for chani = 1:5
        macOut(chani, :) = squeeze(macroDat(chani, :) -...
                                      macroDat(chani+1,:)); 
    end
    macroAvail = true; 
catch
    disp('no macro data')
    macroAvail = false; 
end


    %cleaning spike noise
 switch sessionIDs{sessi}
    case '250811_Dupi_NMH_TPB_1'
        spikeThresh = 10;
        spikeWin = 9; 
        hasEEG = true;
        spikeClean = true; 
    case '250904_OBE_NWU_TI'
        spikeThresh = 50;
        spikeWin = 11; 
        hasEEG = true;
        spikeClean = true; 
    case '250623_DUPI_NMH_KS_2'
        spikeThresh = 15;
        spikeWin = 11; 
        hasEEG = true;
        spikeClean = false; 
    case '250623_Dupi_NMH_KS_1'
        spikeThresh = 15;
        spikeWin = 11; 
        hasEEG = false;
        spikeClean = false; 
    case '250818_Dupi_NMH_JH_2'
        spikeThresh = 20;
        spikeWin = 11; 
        hasEEG = true;
        spikeClean = true; 
    case '250818_Dupi_NMH_JH_1'
        spikeThresh = 20;
        spikeWin = 11; 
        hasEEG = true; 
        spikeClean = true; 
    case '250908_OBE_NWU_AS'
        spikeThresh = 20;
        spikeWin = 11; 
        hasEEG = true;
        spikeClean = false; 
    otherwise
        hasEEG = true; 
        spikeClean = false;
end


%% spike cleaning using prominence detector combined with windowed IC removal
%applied to macro channels only 

if spikeClean & macroAvail
    [b,a] = butter(4, [5,150]/(outDat.fs/2), 'bandpass');
    gammaSig = filtfilt(b,a, macOut')'; 
    [test, prominence] = detect_spikes(macOut,spikeThresh,...
        spikeWin,...
        false, gammaSig); 
        %ICA is on the macro data without bipolar rereference 
        %this allows later rereferencing at will
    out = ica_flag_spikes_targeted(macOut, test, prominence, 'Fs', 500);

   
    outDat.data(end+1:end+5, :,:) = reshape(out.data_clean, 5, T, N);  
    outDat.labels(end+1:end+5) = {'macBP1', 'macBP2', 'macBP3', ...
                                                'macBP4', 'macBP5'};
    outDat.data(end+1, :,:) = reshape(out.mixVector, 1, T, N); 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 
elseif macroAvail
    outDat.data(end+1:end+5, :,:) = reshape(macOut, 5, T, N);
    outDat.labels(end+1:end+5) = {'macBP1', 'macBP2', 'macBP3', ...
                                                'macBP4', 'macBP5'};
    outDat.data(end+1, :,:) = ones(size(outDat.data,[2,3])); 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 


end

%% blink removal using full IC removal across all ephys channels
    if hasEEG
        ephysDat = outDat.data(1:32,:,:);
        [out, badChan, blinkIndicator] = blinkRemoveWrapper(ephysDat,...
                                        outDat.fs);

        outDat.badChans = unique(badChan); 
        outDat.data(1:32,:,:) = out;
        outDat.data(end+1, :,:) = blinkIndicator; 
        outDat.labels{end+1} = "blinkIndicator";
        outDat.blinkRemoval = 1; 
        


       
    end


 

%% pull out respiration data and get breathing summary vars: 
if ~isfield(outDat, 'bmObj')
    %get out the respiration data
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :, :); 
    
    %choose which one looks right
    rspDatz = (rspDat - mean(rspDat, [2,3])) ./ std(rspDat, [], [2,3]); 
    
    idx = rspIDX(sessi); 
    rspDat = squeeze(rspDat(idx, :, :)); 
    rspDatz = squeeze(rspDatz(idx, :,:)); 
    
    %flip signal
    rspDat = rspDat .* rspFlip(sessi);
    
    
    
    
    
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
    
    outDat.bmObj = bmObj; 
    outDat.rspDat = rspDat; 
    
    meanLens = zeros(size(rspDat,2)); 
    for cndi = 1:size(rspDat, 2)
        test = bmObj(bmObj(:,12)==cndi, 2);
        figure
        plot(outDat.tim, rspDat(:, cndi))
        xline(test)
    
        meanLens(cndi) = mean(bmObj(bmObj(:,12)==cndi, 7));
    
    end
    
    %get the target files for the shadow conditions:
    targTraces = zeros(size(rspDat)); 
     for cndi = 3:size(rspDat, 2)
        idx = find(outDat.behDat.order == cndi, 1);
        targFile = outDat.behDat.shadowFile{idx}; 
        if strcmp('NA', targFile)
            targFile = 'audioResp';
        end
        targFile = [sessionIDs{sessi} '_' targFile '_recording.csv']; 
        targFile = ['closed-loop-respiration\data\' ...
                    targFile];
        targFile = readtable([codePre targFile]);
        targFile(targFile.voltage == 0,:) = []; 
        try
            tempo_scale = str2num(outDat.behDat.warp{idx});
        catch
            tempo_scale = outDat.behDat.warp(idx);
        end
    
        targ_len = round( outDat.behDat.trialTim(idx)* ...
                          outDat.behDat.FPS(idx) * 2);
        new_len = round(length(targFile.voltage) / tempo_scale);
        voltages_resampled = interp1(1:length(targFile.voltage),... %og time
                                    targFile.voltage, ... %og signal
                                    1:new_len, 'linear'); %targ time
        loop_start = round(180/tempo_scale); 
        loop_end = length(voltages_resampled) - round(180/tempo_scale);
        loop_segment = voltages_resampled(loop_start:loop_end); 
        loopLen = length(loop_segment); 
        if loopLen > targ_len
            % Truncate the loop_segment to the target length
            voltages = loop_segment(1:targ_len);
        
        elseif loopLen < targ_len
            % Repeat the loop_segment to reach at least the target length
            repeats = ceil(targ_len / length(loop_segment));
            voltages = repmat(loop_segment, 1, repeats);
            voltages = voltages(1:targ_len);
        end
        
        %cut to trialTim length
        timStp = mean(diff(targFile.timestamp));
        tmpTim = timStp:timStp:outDat.behDat.trialTim(idx);
        voltages = voltages(1:length(tmpTim));
        
        %resample to match ephys data: 
        voltages = interp1(tmpTim, ...
                           voltages, ...
                           outDat.tim, 'linear'); 
        
        targTraces(:,cndi) = voltages; 
    
    end
    
    outDat.targTraces = targTraces; 
    outDat.data(end+1, :, :) = targTraces; 
    outDat.labels{end+1} = 'targTrace'; 
    outDat.CSClist{end+1} = 'targTrace'; 
    
    
    % save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
    %                 sessionIDs{sessi} '_breathingPreProc.mat'], ...
    %                 'outDat', "-v7.3")

else
    disp('breathing has already been processed') 

end

   %% get heartbeats: 
        %requires custom heartbeat detection
if ~isfield(outDat, 'RRint')


    %bandpass and z-score
    idx = cellfun(@(x) contains(x, 'ECG'), outDat.labels);
    ECG = outDat.data(idx, :, :); 
    test = reshape(ECG, 3, []); 
    
    d = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', 5, 'HalfPowerFrequency2', 40, ...
    'SampleRate', outDat.fs);
    L = size(outDat.data,2); 
    test = filtfilt(d, test')'; 
    ECG = reshape(test, 3, L, []); 
    
    
    
    %plot for custom algorithm design: 
    ECGz = (ECG - mean(ECG, [2,3])) ./ std(ECG, [], [2,3]); 
    figure; 
    plot(1:L, ECGz(1,:,1), 'color', 'k')
    hold on 
    plot(1:L,ECGz(2,:,1), 'color', 'red')
    plot(1:L,ECGz(3,:,1), 'color', 'green')
    
    
    beatSep = outDat.fs / 10; 
    %participant specific heartbeat analysis:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    heartBeats = cell(size(rspDat, 2),1); 
    switch sessionIDs{sessi}
        case '250818_Dupi_NMH_JH_1'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x>5 & y>4 & z<-.5, ...
                         ECGz(1,3:end,cndi), ECGz(2,1:end-2,cndi), ...
                         ECGz(3,1:end-2,cndi) ) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
        case '250818_Dupi_NMH_JH_2'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x>3 & y<-3 & z>1, ...
                         ECGz(2,1:end-13,cndi), ECGz(3,2:end-12,cndi), ...
                         ECGz(3,14:end,cndi) ) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
        case '250623_DUPI_NMH_KS_2'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z, n) x>.5 & y>1.75 & ...
                                                 z<-2 & n<-1, ...
                         ECGz(1,1:end-9,cndi), ECGz(2,1:end-9,cndi), ...
                         ECGz(3,1:end-9,cndi), ECGz(2,10:end,cndi) ) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
        case '250623_Dupi_NMH_KS_1'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x>2 & y>.75 & ...
                                                 z<-2, ...
                         ECGz(1,1:end-7,cndi), ECGz(2,8:end,cndi), ...
                         ECGz(3,2:end-6,cndi)) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
        case '250908_OBE_NWU_AS'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x>2 & y>5 & z<-4, ...
                         ECGz(1,5:end,cndi), ECGz(2,1:end-4,cndi), ...
                         ECGz(3,2:end-3,cndi) ) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
         case '250723_EEG_NWU_IN'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x<-1 & y>2 & z<-1, ...
                         ECGz(1,5:end,cndi), ECGz(2,3:end-2,cndi), ...
                         ECGz(3,1:end-4,cndi) ) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
         case '250725_EEG_NWU_BN'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x>1 & y<-3 & z>1, ...
                         ECGz(3,12:end,cndi), ECGz(3,1:end-11,cndi), ...
                         ECGz(2,2:end-10,cndi) ) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
         case '250815_EEG_NWU_PP'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y) x<-4 & y>3, ...
                         ECGz(1,1:end,cndi), ECGz(2,1:end,cndi) ) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
         case '250819_EEG_NWU_ZL'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x<-2 & y>4 & z<-2, ...
                         ECGz(1,1:end,cndi), ECGz(2,1:end,cndi), ...
                         ECGz(3,1:end,cndi)) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
         case '250723_EEG_NWU_BK'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x>1.5 & y>2 & z<-3, ...
                         ECGz(1,5:end,cndi), ECGz(2,1:end-4,cndi), ...
                         ECGz(3,3:end-2,cndi)) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
         case '250912_EEG_NWU_JN'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y) x>4 & y<-4, ...
                         ECGz(2,1:end,cndi), ECGz(3,1:end,cndi)) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
         case '250904_OBE_NWU_TI'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x<-2 & y>3 & z<0, ...
                         ECGz(1,1:end-2,cndi), ECGz(2,2:end-1,cndi),...
                         ECGz(3,3:end, cndi)) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
         case '250811_Dupi_NMH_TPB_1'
            for cndi = 1:size(rspDat, 2)
                %get within beat times: 
                test = find(arrayfun(@(x,y,z) x>1 & y<-1 & z<-1, ...
                         ECGz(2,3:end-10,cndi), ECGz(3,1:end-12,cndi),...
                         ECGz(2,13:end, cndi)) ) ;
                test = test(diff(test) > beatSep);
                heartBeats{cndi} = test; 
            end
    
        otherwise
            disp('This participant does not have a set heartbeat alg')
            ekalsjd 
    
    end
    
    outDat.heartBeats = heartBeats; 
    %% with beats detected, establish RR interval timeseries
    
    RRint = zeros(size(rspDat)); 
    for cndi = 1:size(rspDat, 2)
        tim = outDat.tim; 
        beats = heartBeats{cndi}; 
        beatTims = tim(beats); 
        beatDiffs = diff(beatTims);
        beatTims(end) = [];
       
    %  figure; 
    % plot(1:L, ECGz(1,:,cndi), 'color', 'k')
    % hold on 
    % plot(1:L,ECGz(2,:,cndi), 'color', 'red')
    % plot(1:L,ECGz(3,:,cndi), 'color', 'green')
        %account for accidental misses and double counts and interpolate
        %across
        breakVals = [0:.05:3]; 
        counts = arrayfun(@(x,y) sum(beatDiffs>x & beatDiffs<y), ...
            breakVals(1:end-1), breakVals(2:end));
        %locate the mode and find zeros around it to define central dist.
        [~, idx] = max(counts); 
        minVal = breakVals(idx - find(flip(counts(idx-10:idx))==0, 1) + 1); 
        if isempty(minVal)
            minVal = .6; 
        end
        maxVal = breakVals(idx + find(counts(idx:idx+10)==0, 1) - 1); 
    
        %remove bad vals: 
        beatTims(beatDiffs < minVal) = []; 
        beatDiffs(beatDiffs < minVal) = []; 
    
        beatTims(beatDiffs > maxVal) = []; 
        beatDiffs(beatDiffs > maxVal) = []; 
        figure
        histogram(beatDiffs)
         
        RRint(:,cndi) = interp1(beatTims,beatDiffs, tim, 'linear');
    end
    
    figure
    plot(rspDat(:,2))
    yyaxis right
    plot(RRint(:,2))
    
    outDat.data(end+1, :, :) = RRint; 
    outDat.labels{end+1} = 'RRint'; 
    outDat.CSClist{end+1} = 'RRint'; 
    
    % save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
    %                 sessionIDs{sessi} '_breathingPreProc.mat'], ...
    %                 'outDat', "-v7.3")

else
    disp('heart rate has already been processed') 

end


 %% trial epoch all data, eliminating breaths that don't conform properly

 if ~isfield(outDat, 'trialDat')
    %grab 2 s buffer plus breath
    %max of 18 s breath
    L20 = 20*outDat.fs; 
    trialDat = zeros(size(outDat.data,1), ... % channels
                        L20, ...     % 20 seconds
                        size(outDat.bmObj,1)); % breaths
    
    %Good breaths indicator variable for whether a breath is well-behaved: 
    %good fit up to peak from start
    %only one peak and one minimum prior to next breath onset
   
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
%col13: Good Breaths indicator 1 = good 0 = bad
%col14: index
%col15: RRmin
%col16: RRmax
%col17: RRmax - RRmin
    for bb = 1:size(outDat.bmObj, 1)
       bb
        idx = find(outDat.bmObj(bb, 2)<=outDat.tim,1);
        cndi = outDat.bmObj(bb, 12); 
        bL = outDat.bmObj(bb, 7) ; 
        %check breath isn't too near start or end of recording
        %also not more than 18 seconds long
        endPnt = idx + L20  - outDat.fs*2;
        startPnt = endPnt - L20 + 1;  
        if startPnt > 0 && ...
           endPnt < length(outDat.tim) && ...
           bL < 18

          
           curRsp = outDat.rspDat(startPnt:endPnt, cndi); 

           %quality check: 
           %end point should be higher than minimum
           %start point should be higher than minimum
           %one max
           %one min
           smoothRsp = smoothdata(curRsp, 'gaussian', outDat.fs/4);
           bonset = round(outDat.fs*2); 
           boffset = round(outDat.fs*2+bL*outDat.fs); 
           startVal = smoothRsp(bonset); 
           endVal = smoothRsp(boffset);
            
           lowThresh = prctile(smoothRsp(bonset:boffset), 5); 
           highThresh = prctile(smoothRsp(bonset:boffset), 95); 
           
           upPeakIdx = find(smoothRsp(bonset:boffset)<highThresh & ...
                 smoothRsp(bonset+1:boffset+1)>highThresh) ;
            downPeakIdx = find(smoothRsp(bonset:boffset)>highThresh & ...
                               smoothRsp(bonset+1:boffset+1)<highThresh) ;
            
            upTroughIdx = find(smoothRsp(bonset:boffset)<lowThresh & ...
                               smoothRsp(bonset+1:boffset+1)>lowThresh) ;
            downTroughIdx = find(smoothRsp(bonset:boffset)>lowThresh & ...
                                 smoothRsp(bonset+1:boffset+1)<lowThresh) ;
            
            % --- Helper function: true if any pair is 1s apart ---
            tooFar = @(idx) numel(idx) > 1 && any(diff(idx) > outDat.fs);
            
            % --- Reject condition ---
            if startVal < highThresh && ...
               endVal   < highThresh && ...  
               startVal > lowThresh  && ...
               endVal   > lowThresh  && ...
               ~tooFar(upPeakIdx)    && ...
               ~tooFar(downPeakIdx)  && ...
               ~tooFar(upTroughIdx)  && ...
               ~tooFar(downTroughIdx)

               %good breath! 
               trialDat(:,:,bb) = outDat.data(:, startPnt:endPnt, cndi);
               outDat.bmObj(bb, 13) = 1; 

               %get RR variability: 
               outDat.bmObj(bb, 15) = max(trialDat(end, ...
                                            bonset:boffset, bb));
               outDat.bmObj(bb, 16) = min(trialDat(end, ...
                                            bonset:boffset, bb));
               outDat.bmObj(bb, 17) = outDat.bmObj(bb, 15)  -...
                                            outDat.bmObj(bb, 16);

           else
                trialDat(:,:,bb) = NaN; 
                figure
                plot(smoothRsp(bonset:boffset))
                yline([highThresh, lowThresh])
                title(['cndi: ' num2str(cndi) 'bb: ' num2str(bb)])
           end
        end

    end

    outDat.trialDat = trialDat; 
    outDat = rmfield(outDat, "data");
    bmObj = outDat.bmObj; 
    

    % Define column headers
    colNames = { ...
        'onsetY', ...
        'onsetTime', ...
        'peakY', ...
        'peakTime', ...
        'endY', ...
        'endTime', ...
        'length', ...
        'amp', ...
        'peakIdx', ...
        'exhalePeakY', ...
        'exhalePeakTime', ...
        'condition', ...
        'goodBreath', ...
        'index', ...
        'RRmin', ...
        'RRmax', ...
        'RRdif'};
    
    % Convert matrix to table
    bmTable = array2table(bmObj, 'VariableNames', colNames);
    
    % Save table to CSV
    writetable(bmTable, [codePre 'closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '_processedBreathing.csv']);


    save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_breathingPreProc.mat'], ...
                    'outDat', "-v7.3")

else
    disp('all data already epoched') 

end


%% cleaning? 

%7 is length
starts = ones(size(outDat.bmObj,1),1) .* outDat.fs*2; 
stops = round(outDat.bmObj(:,7)*outDat.fs + starts); 
idx = find(cellfun(@(x) contains(x, 'ECG'), outDat.labels));
idx = 1:max(idx); 
[badIDX] = covMatClean_breathing(outDat.trialDat(idx,:,:), starts, stops);


% [data.nbChanOrig, data.nbChanFinal, data.nbTrialOrig, data.nbTrialFinal,behDat, EEG]...
%          = removeNoiseChansVoltAB(EEG, behDat);
% 
% 
% 
% trodeLocs = readtable([codePre ...
%                     'ZelanoLabScripts\myEEGcoords_thetaPhi.csv']);






trialDat = outDat.trialDat; 
bmObj = outDat.bmObj; 

trialDat(:,:,badIDX ) = []; 
bmObj(badIDX , :) = [];
trialDat(:,:,bmObj(:,13) == 0 ) = []; 
bmObj(bmObj(:,13) ==0, :) = [];



idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
rspDat = trialDat(idx, :, :); 


idx = rspIDX(sessi); 
rspDat = squeeze(rspDat(idx, :, :)); 

%flip signal
rspDat = rspDat .* rspFlip(sessi);

figure
[~, order] = sort(bmObj(:,7));
imagesc(.002:.002:20, [], squeeze(rspDat(:,order))')
caxis([-100,100])


%% end

end




%cleaning blinks
        
        % ephysDat = reshape(ephysDat, 32, []);
        % figure; plot(ephysDat([1], :))
        % hold on 
        % plot(ephysDat(32,:))
        % legend({'chan 1', 'chan 32'})
        % blinkChan = input(sprintf(...
        %         'Enter the index of blinkChan (1 or %d), or [] to skip: '...
        %                                         ,size(ephysDat,1)));
        % 
        % 
        % 
        % figure; 
        % imagesc(ephysDat) %check for bad channels overall 
        % caxis([-200, 200])
        % 
        % badChan = input(sprintf(...
        % 'Enter the index of badChans (1..%d), or [] to skip: '...
        %                                         ,size(ephysDat,1)));
        % 
        % figure; 
        % imagesc(ephysDat) %check for bad channels overall 
        % caxis([-1, 1])
        % 
        % badChan = [badChan input(sprintf(...
        % 'Enter the index of badChans (1..%d), or [] to skip: '...
        %                                         ,size(ephysDat,1)))];
        % chanidx = 1:size(ephysDat,1);
        % trainDat = ephysDat; 
        % trainDat(badChan, :) = []; 
        % chanidx(badChan) = []; 
        % 
        % 
        % eyeBlinkDat = ephysDat(blinkChan,:) - mean(trainDat,1); 
        % [b,a] = butter(4, [3,10]/(outDat.fs/2), 'bandpass');
        % eyeBlinkDat = filtfilt(b,a, eyeBlinkDat); 
        % eyeBlinkDat = (eyeBlinkDat - mean(eyeBlinkDat,2)) ./ ...
        %                     std(eyeBlinkDat,[],2);
        % 
        % %check that eyeblink dat looks as expected
        % figure; 
        % plot(eyeBlinkDat')
        % 
        % test = eyeBlinkDat>2;
        % outDat.data(end+1, :,:) = reshape(test, 1, T, N); 
        % outDat.labels{end+1} = "blinkIndicator";
        % 
        % 
        % startIdx = [1:500:length(test)-100000]; 
        % blinkCounts = arrayfun(@(x) sum(test(x:x+100000)), startIdx); 
        % idx = find(blinkCounts>median(blinkCounts)-5 & ...
        %     blinkCounts<median(blinkCounts)+5);
        % 
        % 
        % startIdx = startIdx(idx(1)); 
        % 
        % 
        % 
        % 
        % trainDat = trainDat(:,startIdx:startIdx+100000);
        % 
        % out = ica_blinks(trainDat, 'blinkChan', ...
        %     find(chanidx == blinkChan));
        % 
        % Sclean = out.W * ephysDat(chanidx,:); 
        % Sclean(out.badICs,:) = 0; % removal of blink IC entirely 
        % data_clean = out.A * Sclean;                     % back to channel space
        % newEphys = ephysDat; 
        % newEphys(chanidx,:) = data_clean; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NEW VERSION: 


% % % targTraces = zeros(size(rspDat));  % same size as ephys-aligned matrix
% % % 
% % % for cndi = 3:size(rspDat, 2)
% % % 
% % %     idx = find(outDat.behDat.order == cndi, 1, 'first');
% % %     targFile = outDat.behDat.shadowFile{idx};
% % %     if strcmp(targFile, 'NA'), targFile = 'audioResp'; end
% % % 
% % %     fname = [sessionIDs{sessi} '_' targFile '_recording.csv'];
% % %     fpath = fullfile(codePre, 'closed-loop-respiration', 'data', fname);
% % % 
% % %     T = readtable(fpath);
% % %     % Expect columns named 'timestamp' and 'voltage'
% % %     if ~all(ismember({'timestamp','voltage'}, T.Properties.VariableNames))
% % %         error('Expected columns {timestamp, voltage} in %s', fpath);
% % %     end
% % % 
% % %     % Remove obvious zero fillers, keep columns as double column vectors
% % %     keep = T.voltage ~= 0;
% % %     t_raw = double(T.timestamp(keep));
% % %     v_raw = double(T.voltage(keep));
% % % 
% % %     % tempo scale
% % %     try
% % %         tempo_scale = str2num(outDat.behDat.warp{idx});
% % %     catch
% % %         tempo_scale = outDat.behDat.warp(idx);
% % %     end
% % % 
% % %     % target length in samples, mirroring Python (trialTime * FPS * 2)
% % %     FPS = outDat.behDat.FPS(idx);
% % %     trialTime = outDat.behDat.trialTim(idx);
% % %     targ_len = round(trialTime * FPS * 2);
% % % 
% % %     % --- time-warp (resample) ---
% % %     new_len = max(10, round(numel(v_raw) / tempo_scale));  % guard tiny
% % %     v_resamp = time_warp_resample(t_raw, v_raw, tempo_scale, new_len);
% % % 
% % %     % --- choose loop segment to avoid edges ---
% % %     % emulate Python’s margin ~ 180/tempo_scale (in samples AFTER warp)
% % %     margin = max(1, round(180/tempo_scale));
% % %     if numel(v_resamp) <= 2*margin
% % %         % too short to take margins, just use the whole thing
% % %         loop_seg = v_resamp(:);
% % %     else
% % %         loop_seg = v_resamp( (margin+1) : (numel(v_resamp)-margin) );
% % %     end
% % % 
% % %     % --- tile or truncate to target length ---
% % %     voltages = tile_or_truncate(loop_seg, targ_len);
% % % 
% % %     % --- cut to trialTime length (first half) and align to ephys time grid ---
% % %     % Build the target time vector for the warped “stimulus” at FPS
% % %     t_step = 1/FPS;
% % %     % use exact count equal to half (since you made *2 earlier)
% % %     n_first_half = round(trialTime * FPS);
% % %     voltages = voltages(1:n_first_half);
% % % 
% % %     % resample onto your ephys time base outDat.tim (assume column vector)
% % %     t_stim = (0:n_first_half-1).' * t_step;  % starts at 0
% % %     % ensure outDat.tim is within the stimulus time range
% % %     t_ephys = outDat.tim(:);
% % %     % If outDat.tim extends beyond t_stim, extrapolate or clamp:
% % %     v_aligned = interp1(t_stim, voltages(:), t_ephys, 'linear', 'extrap');
% % % 
% % %     % ensure correct length to assign
% % %     nAssign = min(numel(v_aligned), size(targTraces,1));
% % %     targTraces(1:nAssign, cndi) = v_aligned(1:nAssign);
% % % 
% % % end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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





