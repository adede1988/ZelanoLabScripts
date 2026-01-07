  

clear
datPre  = {fullfile(pwd, 'sample_data')};

sessionIDs = {'exampCueTaskDat'};           
datPrei = [1,1,1,2,2,1,1];       


addpath(genpath(pwd));
%INSERT ADD PATH FOR eeglab2025.0.0 here: 
% addpath(genpath('C:\Users\dtf8829\Documents\eeglab2025.0.0'))

figPath = fullfile(pwd, 'testFigs');
% Load once, reuse
EEGLOC = readtable(fullfile(pwd,'myEEGcoords_thetaPhi.csv'));

%% 
%IGNORE LOOPING STRUCTURE AND FOCUS ON CREATING ONE PIPELINE THAT CAN LATER
%BE WRAPPED IN A FOR LOOP EASILY
sessi = 1; %testing sessi iterator is set here 
% for sessi = 1:length(sessionIDs)

%% data cleaning
%TESTING LOAD SAMPLE DATA: 
 outDat = load(fullfile(pwd, 'sample_data', [sessionIDs{sessi} '.mat'])); 
  

    outDat = downsample_data(outDat, 500);
   


    %cleaning spike noise
 switch sessionIDs{sessi}
    case '250904_OBE_NWU_TI'
        spikeThresh = 50;
        spikeWin = 11; 
        hasEEG = true;
        spikeClean = true; 
    case '230611_OBE_NMH_AZ'
        spikeThresh = 20;
        spikeWin = 5; 
        hasEEG = false;
        spikeClean = true;       
     case '250313_OBE_NMH_CS'
        spikeThresh = 20;
        spikeWin = 5; 
        hasEEG = false;
        spikeClean = false;
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
    case '241017_OBE_NMH_AS'
        spikeThresh = 20;
        spikeWin = 11; 
        hasEEG = false;
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
    case '250310_OBE_NMH_FS'
        spikeThresh = 8;
        spikeWin = 5; 
        hasEEG = false;
        spikeClean = true; 
       

    otherwise
        spikeThresh = 20;
        spikeWin = 11; 
        hasEEG = false;
end

idx = cellfun(@(x) contains(x, 'macro'), outDat.labels);
figure
macroDat = outDat.data(idx, :); 
plot(macroDat(1,:))
hold on 
for ii = 2:6
    plot(macroDat(ii,:)+(ii-1)*50)
end
title('will need to code channel selection if needed')
%% bipolar rereferencing 
macOut = zeros([5, size(macroDat, [2,3])]);
%do bipolar rereferencing 
for chani = 1:5
    macOut(chani, :, :) = squeeze(macroDat(chani, :) -...
                                  macroDat(chani+1,:)); 
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
    chanIDX standardEEGlocs badChan2 out lbls ii chan hasEEG

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
   

    afterWin = 1000; 

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
        case '250313_OBE_NMH_CS'
            thresh = 3000; 
            cuedBackBuff = 150; 
        case '230611_OBE_NMH_AZ' 
            thresh = 3000; 
            cuedBackBuff = 150; 
        case '250310_OBE_NMH_FS' 
            thresh = 3000; 
            cuedBackBuff = 150; 
            afterWin = 2000; 
        case '241017_OBE_NMH_AS' 
            thresh = 3000; 
            cuedBackBuff = 150; 
        case '250811_Dupi_NMH_TPB_1' 
            thresh = 20; 
            cuedBackBuff = 150; 
        case '250904_OBE_NWU_TI' 
            thresh = 5000; 
            cuedBackBuff = 100; 
        case '240923_OBE_NMH_HRM' 
            thresh = 1000; 
            cuedBackBuff = 350; 
        case '250818_Dupi_NMH_JH_1'
            thresh = 4000; 
            cuedBackBuff = 350; 
        case '250818_Dupi_NMH_JH_2'
            thresh = 4000; 
            cuedBackBuff = 350; 
        otherwise
            thresh = 500;
            cuedBackBuff = 150; 
    end

    

    outSniffs = zeros(40, 10); 
    %col 1: sniff onset index
    %col 2: trial number
    %col 3: empty
    %col 4: off from TTL by
    %col 5: trial type (hit/miss/cr/fa
    %col 6: empty 
    %col 7: adjustment for phase align
    oi = 1; %index variable for out sniffs

    used = zeros(size(idx)); 
    for triali = 1:40
        targets = TTLs(triali, 2);
        targets(isnan(targets)) = [];
        targets = sort(targets); 
        sniffi = 1; 
        for tt = 1:length(targets)
            if tt == 1 || tt == length(targets) 
                %first and last sniff! should be after TTL
                startSearch = targets(tt) - 30 - cuedBackBuff; 
                endSearch = targets(tt)+afterWin - 30; %search window of 2s
                val = find(test(startSearch:endSearch-1)<thresh & ... 
                    test(startSearch+1:endSearch) > thresh, 1);
                if tt == 1
                    type = 1; 
                else
                    type = 3; 
                end
             
            else
                %free sniffs! should be before TTL
                startSearch = targets(tt)- 30 - afterWin; %search window of 2s
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
                    outSniffs(oi,4) = - 30 - afterWin +val; 
                else
                    outSniffs(oi,4) = val - 30 - cuedBackBuff;
                end
                %col 5: corect / incorrect behavior 1 / 0
                outSniffs(oi,5) = outDat.behDat.type(triali); 
                %col 6: sniff type 1 = start, 2 = free, 3 = confirm
                outSniffs(oi,6) = type;
                oi = oi + 1; 
            end
        end
    end





    for triali  = 1:40
        figure
        idx = TTLs(triali, 1); 
        endidx = TTLs(triali,3); 
        if isnan(endidx)
            endidx = idx + 7000;
        end
        rspTrial = smoothdata(rspDat(idx:endidx), 'gaussian', 300); 
        testHere = test(idx+30:endidx+30);
     
        plot(rspTrial )
       
        idx = TTLs(triali,2) - TTLs(triali, 1);
        idx(isnan(idx)) = []; 
        xline(idx, 'linewidth', 2)
        title(triali)
        curSniffs = outSniffs(outSniffs(:,2) == triali, 1); 
        if ~isempty(curSniffs)
            xline(curSniffs - TTLs(triali, 1), ...
                'color', 'green', 'linewidth', 2); 
        end
        yyaxis right
        plot(testHere)

    end

     %col 1: sniff onset index
    %col 2: trial number
    %col 3: empty
    %col 4: off from TTL by
    %col 5: trial type (hit/miss/cr/fa
    %col 6: empty 
    %col 7: adjustment for phase align
if size(outSniffs,1) ~= size(outDat.behDat,1)
    error('sniff count does not line up with behavioral data')
end
outDat.moreThan1 = 0; 
outDat.behDat.sniffOnset = outSniffs(:,1); 
outDat.behDat.TTLoffSet = outSniffs(:,4); 



clear cuedBackBuff curSniffs endidx endSearch idx oi outSniffs rspDat ...
    rspTrial sniffi startSearch targets targidx test testHere thresh ...
    triali tt type used val TTLs

%% adjust sniffs and save




 %subjectSpecific Thresholds
switch sessionIDs{sessi}
   
    case '250904_OBE_NWU_TI' 
        adjWin = 300;  
    otherwise
        adjWin = 500; 
end

for sniffi = 1:size(outDat.behDat, 1)
    %loop sniffs
    startIDX = outDat.behDat.sniffOnset(sniffi) - 2*outDat.fs; 
    endIDX = outDat.behDat.sniffOnset(sniffi) + 4*outDat.fs; 
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
    outDat.behDat.adjust(sniffi) = adjust; 
    outDat.behDat.finalOnset(sniffi) = ...
                            outDat.behDat.sniffOnset(sniffi) + adjust; 

   



end

idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
rspDat = outDat.data(idx, :); 

  
    
idx = rspIDX(sessi); 
rspDat = squeeze(rspDat(idx, :)); 
    
%flip signal
rspDat = rspDat .* rspFlip(sessi);



figure; plot(rspDat)
xline(outDat.behDat.finalOnset)


outDat.rspIDX = rspIDX(sessi);
outDat.rspFlip = rspFlip(sessi); 

 if ~exist([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc'], 'dir')
         mkdir([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc']);
 end

 save([datPre{datPrei(sessi)} sessionIDs{sessi} '\preProc\' ...
                    sessionIDs{sessi} '_cueTaskPreProc.mat'], ...
                    'outDat', "-v7.3")

 clear adjust adjWin analytic ans difVals endIDX idx minidx outDat ...
     rspDat rspPhase smthRsp sniffi startIDX test trialDat

% end