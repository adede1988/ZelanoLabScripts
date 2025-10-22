clear

codePre = 'G:\My Drive\GitHub\';
datPre = { 'R:\Neurology\Zelano_Lab\Lab_Common\Dupi\', ... 
           'R:\Neurology\Zelano_Lab\Lab_Common\OBEControl\'};

%prefix index for data folder: 
datPrei = [1,1,1,1,1]; 

sessionIDs = {'250818_Dupi_NMH_JH_1', ... %preprocessed 
               '250623_DUPI_NMH_KS_2',... %preprocessed
               '250623_Dupi_NMH_KS_1',... %preprocessed
                '250818_Dupi_NMH_JH_2',...%preprocessed
                '250811_Dupi_NMH_TPB_1',...%preprocessed
                };   

%there are multiple respiration channels in many recordings
%which one is right for each session: 
rspIDX = [1,1,3,1,1]; 
rspFlip = [1,1,-1,1,1]; %hard code flip

addpath([codePre 'HpcAccConnectivityProject/helperFuncs'])
addpath(genpath([codePre 'myFrequentUse']))
addpath([codePre 'myFrequentUse/export_fig_repo'])
addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))

addpath([codePre 'fieldtrip-20230118'])
addpath([codePre 'emotionDecoding'])
addpath([codePre 'slowBreathing'])

set(0, 'defaultfigurewindowstyle', 'docked')
ft_defaults



for sessi = 1:length(sessionIDs)
%% data load 

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

    idx = cellfun(@(x) contains(x, 'raw_PEAintensityPleasantness'), {datFolders.name});
    idx = find(idx); 
    
    if isscalar(idx)
        %process when two runs need to be stitched together: 

        %load in both files: 
        dat = load([datFolders(idx).folder '/' ...
                     datFolders(idx).name  '/' ...
                     datFolders(idx).name  '.mat']);
        dat = dat.curDat; 
     
        
        %find the behavioral files: 
        behDir = dir([behFold.folder, '/' behFold.name]);
        idx = cellfun(@(x) contains(x, 'PEA'), {behDir.name});
        idx = find(idx);
        if ~isscalar(idx)
            error('wrong number of PEA folders detected')
        end
        behDir = behDir(idx); 
        behDir = dir([behDir.folder '\' behDir.name]);
        idx = cellfun(@(x) contains(x, '.mat'), {behDir.name});
        idx = find(idx); 
        if ~isscalar(idx)
            error('wrong number of behavioral files detected')
        end

        %load the behavioral files
        behDat = [behDir(idx).folder '/' behDir(idx).name];
       
        behDat = load(behDat).outMat;
        behDat = cat(1, behDat{:});
        %fill in empty values
        f = @(x) {735}; 
        emptyIDX = cellfun(@isempty, behDat(:,4));
        behDat(emptyIDX,4) = cellfun(@(x) f(x), behDat(emptyIDX, 4)); 
        emptyIDX = cellfun(@isempty, behDat(:,5));
        behDat(emptyIDX,5) = cellfun(@(x) f(x), behDat(emptyIDX, 5)); 
        behDat = cell2mat(behDat); 
        behDat = array2table(behDat, 'variablenames', ...
                            {'trialNum', 'Odor', 'RT', ...
                             'pleasantness', 'intensity', 'ITI',...
                             'Trial_StartTime', 'Trial_EndTime'});
        
       


        %start making the oudDat struct
        outDat = struct; 
        outDat.labels = dat.outLabs;
        outDat.CSClist = dat.ncslabels; 
        outDat.fs = dat.rawData.fsample; 
        outDat.sessID = sessionIDs{sessi}; 
        outDat.data = dat.rawData.trial{1};
        outDat.behDat = behDat;
        
        %process the TTL pulses for dat: 
        idx = cellfun(@(x) contains(x, 'event'), outDat.labels);
        photoDiode = dat.rawData.trial{1}(idx, :); 

        photoDiode = (photoDiode - mean(photoDiode)) / std(photoDiode);
    

        switch sessionIDs{sessi}
            case '250623_DUPI_NMH_KS_2'
                diodeThresh = -2; 
            otherwise
                diodeThresh = -1.5;
        end

        downs = find(photoDiode(1:end-1) > diodeThresh & ...
            photoDiode(2:end)<diodeThresh);
        ups = find(photoDiode(1:end-1) < diodeThresh & ...
            photoDiode(2:end)>diodeThresh);
        difVals = ups - downs; %difVals is length of TTL pulses
        downs(difVals>1500) = []; 
        ups(difVals>1500) = []; 
        difVals(difVals>1500) = []; 

        downs(difVals<300) = []; 
        ups(difVals<300) = []; 
        difVals(difVals<300) = []; 

        starti = 1; 
        for di = 3:length(downs)
            if downs(di) - downs(di-2) < 3500
                starti = di; 
            end
        end
        downs(1:starti) = []; 
        difVals(1:starti) = []; 

        if length(downs) ~= 45
            error('wrong number of TTLs')
        end
        
        
    
        figure
        plot(photoDiode)
        xline(downs, 'color', 'magenta', 'linewidth', 2)
    

        %store all TTLs into one matrix: 
        %col 1: sniffs
       
        TTLs = nan(45, 1); 
        TTLs(:,1) = downs; 
        TTLs =round(TTLs ./ 4);

   
        
    
    else
        %process for all in two runs: 
        error('what are you supposed to do?')


    end
    
   
   
    %downsample the data
    outDat = downsample_data(outDat, 500);


    clear difVals di downs emptyIDX f idx outMat ...
          photoDiode starti subFolders test ups behDat behDir behFold ...
          datFolders dat

    %% data cleaning

   


     

    %cleaning spike noise
 switch sessionIDs{sessi}
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
    case '250811_Dupi_NMH_TPB_1'
        spikeThresh = 10;
        spikeWin = 9; 
        hasEEG = true;
        spikeClean = true; 
    case '240923_OBE_NMH_HRM'
        spikeThresh = 50;
        spikeWin = 11; 
        hasEEG = false;
        spikeClean = true; 
    otherwise
        spikeThresh = 20;
        spikeWin = 11; 
        hasEEG = false;
 end

  idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
    macroDat = outDat.data(idx, :); 
    macOut = zeros([5, size(macroDat, [2,3])]);
    %do bipolar rereferencing 
    for chani = 1:5
        macOut(chani, :, :) = squeeze(macroDat(chani, :) -...
                                      macroDat(chani+1,:)); 
    end

%% spike cleaning using prominence detector combined with windowed IC removal
%applied to macro channels only 
if spikeClean
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
else
    outDat.data(end+1:end+5, :) = macOut;
    outDat.labels(end+1:end+5) = {'macBP1', 'macBP2', 'macBP3', ...
                                                'macBP4', 'macBP5'};
    outDat.data(end+1, :) = ones(size(outDat.data,2),1); 
    outDat.labels{end+1} = "spikeCleanVec";
    outDat.spikeRemoval = 1; 


end
%% blink removal using full IC removal across all ephys channels
if hasEEG
    %cleaning blinks
    ephysDat = outDat.data(1:32,:); 

    [out, badChan, blinkIndicator] = blinkRemoveWrapper(ephysDat,...
                                    outDat.fs);

    outDat.badChans = badChan; 
    outDat.data(1:32,:) = out; 
    outDat.data(end+1, :,:) = blinkIndicator; 
    outDat.labels{end+1} = "blinkIndicator";
    outDat.blinkRemoval = 1; 
end


 clear a b badChan blinkIndicator chani ephysDat gammaSig hasEEG macOut ...
     macroDat out prominence spikeClean spikeThresh spikeWin
  
%% respiration data handling to identify sniffs

    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx, :); 
    
  
    
    idx = rspIDX(sessi); 
    rspDat = squeeze(rspDat(idx, :)); 
    
    %flip signal
    rspDat = rspDat .* rspFlip(sessi);
    
    
   

    rspTrial = smoothdata(rspDat, 'gaussian', 300); 
    test = (rspTrial(30:end) - rspTrial(1:end-29)).^2 .* ...
    rspTrial(30:end); 
    test(test<0) = 0; 
    test(test>10000) = 10000; 
    test = smoothdata(test, 'gaussian', 500);
   


    %subjectSpecific Thresholds
    switch sessionIDs{sessi}
        case '250623_DUPI_NMH_KS_2' 
            thresh = 3000; 
            cuedBackBuff = 350; 
        case '250623_Dupi_NMH_KS_1' 
            thresh = 5000; 
            cuedBackBuff = 350; 
        case '250908_OBE_NWU_AS' 
            thresh = 3000; 
            cuedBackBuff = 150; 
        case '250818_Dupi_NMH_JH_2'
            thresh = 5000; 
            cuedBackBuff = 350; 
        case '250811_Dupi_NMH_TPB_1' 
            thresh = 3000; 
            cuedBackBuff = 150; 
        case '250904_OBE_NWU_TI' 
            thresh = 5000; 
            cuedBackBuff = 100; 
        case '240923_OBE_NMH_HRM' 
            thresh = 1000; 
            cuedBackBuff = 350; 
        otherwise
            thresh = 500;
            cuedBackBuff = 300; 
    end

    

    outSniffs = zeros(size(outDat.behDat, 1), 10); 
    %col 1: sniff onset index
    %col 2: trial number
    %col 3: empty
    %col 4: off from TTL by
    %col 5: trial type (hit/miss/cr/fa
    %col 6: empty 
    %col 7: adjustment for phase align
    oi = 1; %index variable for out sniffs

    used = zeros(size(idx)); 
    for triali = 1:size(outDat.behDat, 1)
        targets = TTLs(triali, 1);
        targets(isnan(targets)) = [];
        targets = sort(targets); 
        sniffi = 1; 
        for tt = 1:length(targets)
            if tt == 1 || tt == length(targets) 
                %first and last sniff! should be after TTL
                startSearch = targets(tt) - 30 - cuedBackBuff; 
                endSearch = targets(tt)+1000 - 30; %search window of 2s
                val = find(test(startSearch:endSearch-1)<thresh & ... 
                    test(startSearch+1:endSearch) > thresh, 1);
                if tt == 1
                    type = 1; 
                else
                    type = 3; 
                end
             
            else
                %free sniffs! should be before TTL
                startSearch = targets(tt)- 30 - 1000; %search window of 2s
                endSearch = targets(tt) -30 + 300; 

                val = find(test(startSearch:endSearch-1)<thresh & ... 
                    test(startSearch+1:endSearch) > thresh);
                val(val<16) = 16; 
                val(val>1285) = 1285; 
                if length(val)>1
                    testSpace = test(startSearch:endSearch-1);
                    sizes = testSpace(val+15) - testSpace(val-15); 
                    [~, maxSniff] = max(sizes); 
                    val = val(maxSniff); 
                end
                type = 2; 
            end
            
                
                %was a sniff found? if yes then find overall idx
                %val is how far from the start search the sniff happened
            if ~isempty(val)
                targidx = val + startSearch; 
              
                %col 1: sniff onset index
                outSniffs(oi,1) = targidx; 
                %col 2: trial number
                outSniffs(oi,2) = triali; 
                %col 3: sniff within trial number
                outSniffs(oi,3) = sniffi; 
                sniffi = sniffi + 1; 
                %col 4: off from TTL by
                if type == 2
                    outSniffs(oi,4) = - 30 - 1000 +val; 
                else
                    outSniffs(oi,4) = val - 30 - cuedBackBuff;
                end
                %col 5: corect / incorrect behavior 1 / 0
                outSniffs(oi,5) = 0; 
                %col 6: sniff type 1 = start, 2 = free, 3 = confirm
                outSniffs(oi,6) = type;
                oi = oi + 1; 
            end
        end
    end





    for triali  = 1:size(outDat.behDat,1)
        figure
        idx = TTLs(triali, 1) - 1000; 
        endidx = idx + 2000; 
        rspTrial = smoothdata(rspDat(idx:endidx), 'gaussian', 300); 
        testHere = test(idx+30:endidx+30);
     
        plot(rspTrial )
       
        
        xline(1000, 'linewidth', 2)
        title(triali)
        curSniffs = outSniffs(outSniffs(:,2) == triali, 1); 
        if ~isempty(curSniffs)
            xline(curSniffs - TTLs(triali, 1) + 1000, ...
                'color', 'green', 'linewidth', 2); 
        end
        yyaxis right
        plot(testHere)

    end

outDat.sniffInfo = outSniffs; 

clear rspDat rspTrial startSearch targets targidx test testHere ...
    thresh triali tt TTLs type used val curSniffs cuedBackBuff endidx ...
    endSearch idx oi outSniffs
outDat.sniffInfo(outDat.sniffInfo(:,1) == 0, :) = [];
%% epoch the data! 

%channels X 6 seconds of time X sniffs
outDat.tim = -2:1/outDat.fs:4; 
trialDat = zeros(size(outDat.data,1), length(outDat.tim), ...
    size(outDat.sniffInfo,1));


 %subjectSpecific Thresholds
switch sessionIDs{sessi}
   
    case '250904_OBE_NWU_TI' 
        adjWin = 300;  
    otherwise
        adjWin = 500; 
end

for sniffi = 1:size(outDat.sniffInfo, 1)
    %loop sniffs
    startIDX = outDat.sniffInfo(sniffi,1) - 2*outDat.fs; 
    endIDX = outDat.sniffInfo(sniffi,1) + 4*outDat.fs; 
    test = outDat.data(:,startIDX:endIDX); 
    
    %pull trial respiration data: 
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = test(idx, :); 



    idx = rspIDX(sessi); 
    rspDat = squeeze(rspDat(idx, :)); 

    %flip signal
    rspDat = rspDat .* rspFlip(sessi);
    figure
    plot(outDat.tim, rspDat)
    title(sniffi)

    analytic = hilbert(rspDat); 
    analytic = smoothdata(analytic, 1, 'gaussian', 500);
    analytic = lowpass(analytic, 1, outDat.fs); 
    rspPhase = angle(analytic);

    yyaxis right
    plot(outDat.tim, rspPhase)

   %detect most likely peak of the sniff: 
    adjust = find(rspPhase(1:end-1) < 0 & rspPhase(2:end)>0);

    
    if length(adjust) > 1
        adjust = adjust(find(adjust>1000, 1));
    end
    xline(outDat.tim(adjust))
    %go backwards to find the onset
    adjust = adjust - 10;
    smthRsp = smoothdata(rspDat, 'gaussian', 500); 
    difVals = smthRsp(adjust-adjWin:adjust-20) - ...
        smthRsp(adjust-(adjWin-20):adjust);
    yyaxis left
    hold on 
    plot(outDat.tim(adjust-(adjWin-10):adjust-10), difVals*5)

    [~, minidx] = min(difVals); 

    xline(outDat.tim(adjust-(adjWin+20) + minidx))
    


    adjust = round(adjust-520 + minidx - outDat.fs*2);
    test = outDat.data(:,adjust+startIDX:endIDX+adjust); 
    outDat.sniffInfo(sniffi, 7) = adjust; 
    outDat.sniffInfo(sniffi, 8) = outDat.sniffInfo(sniffi, 1) + adjust; 

   

    trialDat(:,:,sniffi) = test;


end

idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
rspDat = trialDat(idx, :, :); 

  
    
idx = rspIDX(sessi); 
rspDat = squeeze(rspDat(idx, :, :)); 
    
%flip signal
rspDat = rspDat .* rspFlip(sessi);

figure; imagesc(outDat.tim, [], rspDat')

figure; plot(rspDat)

outDat.trialDat = trialDat; 
outDat = rmfield(outDat, 'data');

outDat.rspIDX = rspIDX(sessi);
outDat.rspFlip = rspFlip(sessi); 

 if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc'], 'dir')
         mkdir([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc']);
 end

 save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_PEA_threshold_preproc.mat'], ...
                    'outDat', "-v7.3")

end









