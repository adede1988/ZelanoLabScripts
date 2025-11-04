


%custom function to stitch together breathing recordings into one file

clear

codePre = 'G:\My Drive\GitHub\';
datPre = { 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ... 
           'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\',...
           'R:\Neurology\Zelano_Lab\Lab_Common\AllStudyData\EEGbreathing\'};

%prefix index for data folder: 
datPrei = [1,1,1,2,3,3,3,3,3,3,2,1,1]; 

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
                '250811_Dupi_NMH_TPB_1'};

%there are multiple respiration channels in many recordings
%which one is right for each session: 
rspIDX = [3,3,3,3,1,1,1,1,1,1,3,3,3]; 
rspFlip = [-1,-1,-1,-1,-1,-1,-1,1,-1,1,1,1,1]; %hard code flip

% addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
% addpath(genpath([codePre 'myFrequentUse']))
% addpath([codePre 'myFrequentUse/export_fig_repo'])

addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))
% addpath([codePre 'fieldtrip-20230118'])
% addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])
addpath([codePre 'ZelanoLabScripts'])

set(0, 'defaultfigurewindowstyle', 'docked')
% ft_defaults

for sessi = 1:length(sessionIDs)

%% custom import for different participants: 

%check for pre existing processing: 
% if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' 'TESTTESTTESTTESTTEST' ...
%                 sessionIDs{sessi} '_breathingPreProc.mat'], 'file')




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
    
    %put all the data into Chan X time with conditions concatenated
    %together

    data = zeros(size(outDat.data,1), prod(size(outDat.data,[2,3])));
    for ii = 1:size(outDat.data,3)
        data(:,1+(ii-1)*600000:(ii*600000)) = outDat.data(:,:,ii); 
    end
    
    outDat.data = data; 
    % downsample and bandpass: 
    outDat = downsample_data(outDat, 500);
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
    outDat.preProcScript = 'BreathingTaskPreProc.m'; 
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
    % save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
    %                 sessionIDs{sessi} '_breathingPreProc.mat'], ...
    %                 'outDat', "-v7.3")
%if there is pre existing processing, then load it: 
% else
%     outDat =  load([datPre{datPrei(sessi)} sessionIDs{sessi} ...
%                           '\preProc\' sessionIDs{sessi} ...
%                           '_breathingPreProc.mat']).outDat; 
% 
% end

clear behDat dat dat2 data di ii photoDiodeDat tim tmp TTLs secondBreaks

    %% combined flow from here: 

    %at this point there should be an outDat struct with the following: 
    %data: channels X time X conditions
    %tim: 1Xtime vector in seconds
    %behDat: tidy table of behavioral responses to emotion questions
    %labels: 1Xchannels cell array of human readable channel labels
    %CSClist: 1Xchannels cell array of original CSC labels from neuralynx
    %fs: scalar sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% data cleaning
X = outDat.data;  

try
    idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
    macroDat = X(idx, :); 
    idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
    figure
    macroDat = outDat.data(idx, :); 
    plot(macroDat(1,:))
    hold on 
    for ii = 2:6
        plot(macroDat(ii,:)+(ii-1)*50)
    end
    title('will need to code channel selection if needed')
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

%% EEG channel flagging, interpolation, blink removal, laplacian

if hasEEG
    standardEEGlocs = readtable([codePre ...
        'ZelanoLabScripts/myEEGcoords_thetaPhi.csv']);

    %check electrode names
    for chan = 1:32
        if ~strcmp(standardEEGlocs.Label(chan), outDat.labels(chan))
            error('EEG channels not labeled as expected')
        end
    end

    [badTS, badChans] = ...
            removeNoiseChansVolt(outDat.data(1:32,:), outDat.fs);
    chanIDX = 1:32; 
    if ismember(1, badChans)
      chanIDX(badChans(2:end)) = [];
    else
      chanIDX(badChans) = [];
      
    end
    ephysDat = outDat.data(1:32,:);
    [out, badChan2, blinkIndicator] = blinkRemoveWrapper(...
                                                ephysDat(chanIDX,:),...
                                                outDat.fs);

    badChans = [badChans; chanIDX(badChan2)]; 
    ephysDat(chanIDX,:) = out; 
    phi = standardEEGlocs.Phi .*pi ./ 180; 
    theta = standardEEGlocs.Theta .*pi ./ 180; 
    X = cos(phi) .* sin(theta);
    Y = sin(phi) .* sin(theta); 
    Z = cos(theta); 

    ephysDat = interpolate_perrinX(ephysDat,X,Y,Z,badChans);

    ephysDat = ephysDat - mean(ephysDat,1); 

    dataLap = laplacian_perrinX(ephysDat, X, Y, Z); 

    lbls = string(outDat.labels(1:32));
    lbls = reshape(lbls, [length(lbls), 1]);

    outDat.eegLocs = table(lbls, 'VariableNames', {'labels'}); 
    outDat.eegLocs.X = X; 
    outDat.eegLocs.Y = Y; 
    outDat.eegLocs.Z = Z; 
    outDat.eegLocs.theta = theta; 
    outDat.eegLocs.phi = phi; 
    outDat.badChans = outDat.labels(badChans); 
    outDat.dataLap = dataLap; 
    outDat.data(1:32,:) = ephysDat; 
    outDat.data(end+1, :,:) = blinkIndicator; 
    outDat.labels{end+1} = "blinkIndicator";
    outDat.data(end+1, :,:) = badTS; 
    outDat.labels{end+1} = "badTS";
    outDat.EEGInterpolation = 1;
    outDat.EEGCleaning = 1;
    outDat.blinkRemoval = 1; 


end
clear X Y Z theta phi ephysDat dataLap blinkIndicator badTS badChans ...
    chanIDX standardEEGlocs badChan2 out lbls ii chan hasEEG chani idx

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

   
    outDat.data(end+1:end+5, :) = out.data_clean; 
    outDat.labels(end+1:end+5) = {'macBP1', 'macBP2', 'macBP3', ...
                                                'macBP4', 'macBP5'};
    outDat.data(end+1, :) = out.mixVector; 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 
elseif macroAvail
    outDat.data(end+1:end+5, :) = macOut;
    outDat.labels(end+1:end+5) = {'macBP1', 'macBP2', 'macBP3', ...
                                                'macBP4', 'macBP5'};
    outDat.data(end+1, :) = ones(size(outDat.data,[2]),1); 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 


end


clear a b chani gammaSig idx macOut macroAvail macroDat out prominence ...
    spikeClean spikeThresh spikeWin test

 

%% pull out respiration data and get breathing summary vars: 
if ~isfield(outDat, 'bmObj')
    %get out the respiration data
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :); 
    
    %choose which one looks right
    % rspDatz = (rspDat - mean(rspDat, [2,3])) ./ std(rspDat, [], [2,3]); 
    
    idx = rspIDX(sessi); 
    rspDat = squeeze(rspDat(idx, :)); 
    % rspDatz = squeeze(rspDatz(idx, :,:)); 
    
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
    
    bmObj = breathTemplates4(rspDat, outDat.fs);
    bmObj(:, 12) = 1; 
    TaskBreaks = 0:300:max(outDat.behDat.order)*300-10; 
    for cndi = 2:length(TaskBreaks)
       
        bmObj(bmObj(:,2)>TaskBreaks(cndi),12) = cndi; 
    
    end
    outDat.moreThan1 = 1; 
    tmpBehDat = outDat.behDat; 

    outDat.behDat = table(); 
    %col 2: onset tim
    tim = (1:size(outDat.data,2)) / outDat.fs; 
    idx = arrayfun(@(x) find(x<=tim, 1), bmObj(:,2)); 
    outDat.behDat.sniffOnset = idx; 
    outDat.behDat.finalOnset = idx; 
    %col12: condition
    outDat.behDat.condition = bmObj(:,12); 
    %col 1: onset Y value
    outDat.behDat.Yonset = bmObj(:,1); 
    %col 3: peak Y value
    outDat.behDat.inhaleMax = bmObj(:,3); 
    %col 4: peak tim
    idx = arrayfun(@(x) find(x<=tim, 1), bmObj(:,4)); 
    outDat.behDat.inMaxTim = idx; 
    %col 5: end Y value
    outDat.behDat.Yend = bmObj(:,5); 
    %col 6: end tim
    idx = arrayfun(@(x) find(x<=tim, 1), bmObj(:,6)); 
    outDat.behDat.endTim = idx; 
    %col 7: length (end tim - onset tim)
    outDat.behDat.length = bmObj(:,7); 
    %col 8: amp (peak Y - avg of two ends)
    outDat.behDat.amp = bmObj(:,8); 
    %col10: exhale peak Y value
    outDat.behDat.exhaleMin = bmObj(:,10); 
    %col11: exhale peak tim
    outDat.behDat.exMinTim = bmObj(:,11); 
    %col14: index
    outDat.behDat.index = bmObj(:,14); 
    
    %integrate emotion data into the respiration data in general 
    Qs = unique(tmpBehDat.Q_short); 
    
    for cndi = 1:max(outDat.behDat.condition)
        idx = find(outDat.behDat.condition == cndi); 
        tmp = tmpBehDat(tmpBehDat.order == cndi,:);
        outDat.behDat.task(idx) = tmp.task(1); 
        outDat.behDat.noseMouth(idx) = tmp.noseMouth(1);
        outDat.behDat.shadowFile(idx) = tmp.shadowFile(1);
        outDat.behDat.warp(idx) = tmp.warp(1);
        for q = 1:length(Qs)
            ii = find(cellfun(@(x) strcmp(Qs{q}, x), tmp.Q_short));
            outDat.behDat.([Qs{q} '_' tmp.type{ii}])(idx) = tmp.rsp(ii); 
        end
    end

    cndi = 0;
    baseEmotion = table; 
    tmp = tmpBehDat(tmpBehDat.order == cndi,:);
    baseEmotion.task = tmp.task(1); 
    baseEmotion.noseMouth = tmp.noseMouth(1);
    baseEmotion.shadowFile = tmp.shadowFile(1);
    baseEmotion.warp = tmp.warp(1);
    for q = 1:length(Qs)
        ii = find(cellfun(@(x) strcmp(Qs{q}, x), tmp.Q_short));
        baseEmotion.([Qs{q} '_' tmp.type{ii}]) = tmp.rsp(ii); 
    end
    
    
    outDat.baseEmotion = baseEmotion; 
    outDat.rspIdx = rspIDX(sessi);
    outDat.rspFlip = rspFlip(sessi); 
    
    % meanLens = zeros(size(rspDat,2)); 
    % for cndi = 1:size(rspDat, 2)
    %     test = bmObj(bmObj(:,12)==cndi, 2);
    %     figure
    %     plot(outDat.tim, rspDat(:, cndi))
    %     xline(test)
    % 
    %     meanLens(cndi) = mean(bmObj(bmObj(:,12)==cndi, 7));
    % 
    % end
    
    %get the target files for the shadow conditions:
    targTraces = zeros(length(TaskBreaks), outDat.fs * 300); 
     for cndi = 3:max(outDat.behDat.condition)
        idx = find(tmpBehDat.order == cndi, 1);
        targFile = tmpBehDat.shadowFile{idx}; 
        if strcmp('NA', targFile)
            targFile = 'audioResp';
        end
        targFile = [sessionIDs{sessi} '_' targFile '_recording.csv']; 
        targFile = ['closed-loop-respiration\data\' ...
                    targFile];
        targFile = readtable([codePre targFile]);
        targFile(targFile.voltage == 0,:) = []; 
        try
            tempo_scale = str2num(tmpBehDat.warp{idx});
        catch
            tempo_scale = tmpBehDat.warp(idx);
        end
    
        targ_len = round( tmpBehDat.trialTim(idx)* ...
                          tmpBehDat.FPS(idx) * 2);
        new_len = round(length(targFile.voltage) / tempo_scale);

        L = length(targFile.voltage); 
        timRec = 1/L:1/L:1;
        timGoal= 1/new_len:1/new_len:1; 
        voltages_resampled = interp1(timRec,... %og time
                                    targFile.voltage, ... %og signal
                                    timGoal, 'linear'); %targ time
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
        tmpTim = timStp:timStp:tmpBehDat.trialTim(idx);
        voltages = voltages(1:length(tmpTim));
        
        %resample to match ephys data: 
        voltages = interp1(tmpTim, ...
                           voltages, ...
                           1/outDat.fs:1/outDat.fs:300, 'linear');  %FLAG HERE!!!!!!!
        
        targTraces(cndi,:) = voltages; 
    
    end
    targTraces = targTraces'; 
    outDat.data(end+1, :) = targTraces(:); 
    outDat.labels{end+1} = 'targTrace'; 
    
    figure
    plot(rspDat)
    yyaxis right
    plot(targTraces(:))
    
    % save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
    %                 sessionIDs{sessi} '_breathingPreProc.mat'], ...
    %                 'outDat', "-v7.3")

else
    disp('breathing has already been processed') 

end

clear bmObj baseEmotion cndi idx ii L loop_end loop_segment loop_start ...
    loopLen new_len photoDiode q Qs repeats rspDat targ_len targFile ...
    targTraces TaskBreaks tempo_scale tim timGoal timRec timStp tmp ...
    tmpBehDat tmpTim voltages voltages_resampled

   %% get heartbeats: 
        %requires custom heartbeat detection
if ~isfield(outDat, 'RRint')


    %bandpass and z-score
    idx = cellfun(@(x) contains(x, 'ECG'), outDat.labels);
    ECG = outDat.data(idx, :); 
    
    d = designfilt('bandpassiir', 'FilterOrder', 4, ...
    'HalfPowerFrequency1', 5, 'HalfPowerFrequency2', 40, ...
    'SampleRate', outDat.fs);
  
    ECG = filtfilt(d, ECG')'; 
    
    
    
    %plot for custom algorithm design: 
    ECGz = (ECG - mean(ECG, 2)) ./ std(ECG, [], 2); 
    figure; 
    plot(ECGz(1,:), 'color', 'k')
    hold on 
    plot(ECGz(2,:), 'color', 'red')
    plot(ECGz(3,:), 'color', 'green')
    
    
    beatSep = outDat.fs / 20; 
    %participant specific heartbeat analysis:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch sessionIDs{sessi}
        case '250818_Dupi_NMH_JH_1'
            
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x>5 & y>4 & z<-.5, ...
                     ECGz(1,3:end), ECGz(2,1:end-2), ...
                     ECGz(3,1:end-2) ) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         
        case '250818_Dupi_NMH_JH_2'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x>3 & y<-3 & z>1, ...
                     ECGz(2,1:end-13), ECGz(3,2:end-12), ...
                     ECGz(3,14:end) ) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 

        case '250623_DUPI_NMH_KS_2'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z, n) x>.5 & y>1.75 & ...
                                             z<-2 & n<-1, ...
                     ECGz(1,1:end-9), ECGz(2,1:end-9), ...
                     ECGz(3,1:end-9), ECGz(2,10:end) ) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
        case '250623_Dupi_NMH_KS_1'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x>2 & y>.75 & ...
                                             z<-2, ...
                     ECGz(1,1:end-7), ECGz(2,8:end), ...
                     ECGz(3,2:end-6)) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
        case '250908_OBE_NWU_AS'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x>2 & y>5 & z<-4, ...
                     ECGz(1,5:end), ECGz(2,1:end-4), ...
                     ECGz(3,2:end-3) ) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         case '250723_EEG_NWU_IN'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x<-1 & y>2 & z<-1, ...
                     ECGz(1,5:end), ECGz(2,3:end-2), ...
                     ECGz(3,1:end-4) ) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         case '250725_EEG_NWU_BN'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x>1 & y<-3 & z>1, ...
                     ECGz(3,12:end), ECGz(3,1:end-11), ...
                     ECGz(2,2:end-10) ) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         case '250815_EEG_NWU_PP'
            %get within beat times: 
            test = find(arrayfun(@(x,y) x<-4 & y>3, ...
                     ECGz(1,1:end), ECGz(2,1:end) ) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         case '250819_EEG_NWU_ZL'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x<-2 & y>4 & z<-2, ...
                     ECGz(1,1:end), ECGz(2,1:end), ...
                     ECGz(3,1:end)) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         case '250723_EEG_NWU_BK'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x>1.5 & y>2 & z<-3, ...
                     ECGz(1,5:end), ECGz(2,1:end-4), ...
                     ECGz(3,3:end-2)) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         case '250912_EEG_NWU_JN'
            %get within beat times: 
            test = find(arrayfun(@(x,y) x>4 & y<-4, ...
                     ECGz(2,1:end), ECGz(3,1:end)) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         case '250904_OBE_NWU_TI'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x<-2 & y>3 & z<0, ...
                     ECGz(1,1:end-2), ECGz(2,2:end-1),...
                     ECGz(3,3:end, cndi)) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
         case '250811_Dupi_NMH_TPB_1'
            %get within beat times: 
            test = find(arrayfun(@(x,y,z) x>1 & y<-1 & z<-1, ...
                     ECGz(2,3:end-10), ECGz(3,1:end-12),...
                     ECGz(2,13:end, cndi)) ) ;
            test = test(diff(test) > beatSep);
            heartBeats = test; 
    
        otherwise
            disp('This participant does not have a set heartbeat alg')
            ekalsjd 
    
    end
    
    outDat.heartBeats = heartBeats; 

    clear beatSep d ECG ECGz idx heartBeats test

    %% with beats detected, establish RR interval timeseries
    
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :); 
    
    %choose which one looks right
    % rspDatz = (rspDat - mean(rspDat, [2,3])) ./ std(rspDat, [], [2,3]); 
    
    idx = rspIDX(sessi); 
    rspDat = squeeze(rspDat(idx, :)); 
    % rspDatz = squeeze(rspDatz(idx, :,:)); 
    
    %flip signal
    rspDat = rspDat .* rspFlip(sessi);
    tim = 1/outDat.fs:1/outDat.fs:size(outDat.data,2)/outDat.fs; 
    beatTims = tim(outDat.heartBeats); 
    beatDiffs = diff(beatTims);
    beatTims(end) = [];
       
   
    %account for accidental misses and double counts and interpolate
    %across
    breakVals = [0:.05:10]; 
    counts = arrayfun(@(x,y) sum(beatDiffs>x & beatDiffs<y), ...
    breakVals(1:end-1), breakVals(2:end));
    %locate the mode and find zeros around it to define central dist.
    [~, idx] = max(counts); 
    minVal = breakVals(idx - find(flip(counts(1:idx))<5, 1) + 1); 
    if isempty(minVal)
        minVal = .6; 
    end
    maxVal = breakVals(idx + find(counts(idx:end)<5, 1) - 1); 

    %remove bad vals: 
    beatTims(beatDiffs < minVal) = []; 
    beatDiffs(beatDiffs < minVal) = []; 

    beatTims(beatDiffs > maxVal) = []; 
    beatDiffs(beatDiffs > maxVal) = []; 
    figure
    histogram(beatDiffs)
         
    RRint = interp1(beatTims,beatDiffs, tim, 'linear');
    
    
    figure
    plot(rspDat)
    yyaxis right
    plot(RRint)
    
    outDat.data(end+1, :) = RRint; 
    outDat.labels{end+1} = 'RRint'; 

    
    % save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
    %                 sessionIDs{sessi} '_breathingPreProc.mat'], ...
    %                 'outDat', "-v7.3")

else
    disp('heart rate has already been processed') 

end

clear beatDiffs beatTims breakVals counts idx maxVal minVal RRint ...
    rspDat tim

 %% trial epoch all data, eliminating breaths that don't conform properly

 if ~isfield(outDat, 'trialDat')
    %grab 2 s buffer plus breath
    %max of 18 s breath
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :); 
    
    %choose which one looks right
    % rspDatz = (rspDat - mean(rspDat, [2,3])) ./ std(rspDat, [], [2,3]); 
    
    idx = rspIDX(sessi); 
    rspDat = squeeze(rspDat(idx, :)); 
    % rspDatz = squeeze(rspDatz(idx, :,:)); 
    
    %flip signal
    rspDat = rspDat .* rspFlip(sessi);
    idx = cellfun(@(x) contains(x, 'RRint'), outDat.labels);
    rrDat = outDat.data(idx, :); 
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
    for bb = 1:size(outDat.behDat, 1)
       bb
        idx = outDat.behDat.finalOnset(bb);
        bL = outDat.behDat.length(bb) ; 
        %check breath isn't too near start or end of recording
        %also not more than 18 seconds long
        endPnt = idx + outDat.fs*20  - outDat.fs*2;
        startPnt = endPnt - outDat.fs*20 + 1;  
        if startPnt > 0 && ...
           endPnt < length(outDat.tim) && ...
           bL < 18

          
            curRsp = rspDat(startPnt:endPnt); 

           %quality check: 
           %end point should be higher than minimum
           %start point should be higher than minimum
           %one max
           %one min
           smoothRsp = smoothdata(curRsp, 'gaussian', outDat.fs/2);
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
               outDat.behDat.goodBreath(bb) = 1; 

               %get RR variability: 
                %col15: RRmin
                %col16: RRmax
                %col17: RRmax - RRmin
               outDat.behDat.maxRR(bb) = max(rrDat( ...
                                            bonset:boffset));
               outDat.behDat.minRR(bb) = min(rrDat( ...
                                            bonset:boffset));
               outDat.behDat.RR_max_min(bb) = outDat.behDat.maxRR(bb) -...
                                            outDat.behDat.minRR(bb);

           else
                outDat.behDat.goodBreath(bb) = 0; 

               %get RR variability: 
                %col15: RRmin
                %col16: RRmax
                %col17: RRmax - RRmin
               outDat.behDat.maxRR(bb) = max(rrDat( ...
                                            bonset:boffset));
               outDat.behDat.minRR(bb) = min(rrDat( ...
                                            bonset:boffset));
               outDat.behDat.RR_max_min(bb) = outDat.behDat.maxRR(bb) -...
                                            outDat.behDat.minRR(bb);
                figure
                plot(smoothRsp(bonset:boffset))
                yline([highThresh, lowThresh])
                title(['bb: ' num2str(bb)])
           end
        end

    end

   
    
    % Save table to CSV
    writetable(outDat.behDat, [codePre 'closed-loop-respiration\processedBehavior\' ...
                    sessionIDs{sessi} '_processedBreathing.csv']);


    save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_breathingPreProc.mat'], ...
                    'outDat', "-v7.3")

else
    disp('all data already epoched') 

end


idx = outDat.behDat.sniffOnset(outDat.behDat.goodBreath == 1); 
test = arrayfun(@(x) rspDat(x-1000:x+5000), idx, 'uniformoutput', false);
test = cat(1, test{:});

lenVals = outDat.behDat.length(outDat.behDat.goodBreath == 1); 
figure
[~, order] = sort(lenVals);
imagesc(-2:.002:10, [], test(order,:))
caxis([-100,100])

test = arrayfun(@(x) rrDat(x-1000:x+5000), idx, 'uniformoutput', false);
test = cat(1, test{:});

test = (test - mean(test, 2)) ./ std(test, [], 2);


figure
[~, order] = sort(lenVals);
imagesc(-2:.002:10, [], test(order,:))


clear bb bL boffset bonset curRsp downPeakIdx downTroughIdx endPnt endVal...
    highThresh idx lowThresh lenVals order rrDat rspDat smoothRsp ...
    startPnt startVal test tooFar upPeakIdx upTroughIdx


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





