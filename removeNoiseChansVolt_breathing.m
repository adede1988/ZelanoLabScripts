function [nbChanOrig, nbChanFinal, nbTrialOrig, nbTrialFinal, behDat, EEG] = ...
    removeNoiseChansVolt_breathing(EEG, behDat)

%input: 
%       EEG             EEGlab struct with one subject's data
%       path            file path of where the data came from

%output: 
%       nbChanOrig      original number of channels in the dataset
%       nbChanFinal     final number of channels after cleaning
%       nbTrialOrig     original number of trials in the dataset
%       nbTrialFinal    final number of trials in the dataset

%this function calculates the maximum observed deflection in a sliding 80ms
%window across each trial for each channel in the data
%noise channels are defined as channels that have over 50% of their trials
%with max deflections of greater than 100uV
%noise trials are defined as trials that have over 25% of their channels
%with max deflections of greater than 100uV deflections

%in general, noise channels are removed first. However, if over 25% of
%channels are flagged as noise, then trials with 75% of channels with max
%deflections of greater than 100uV are removed first.

%written by Adam Dede (adam.osman.dede@gmail.com)
%summer 2022
% EEG.zDat = EEG.data; 
% for chan = 1:size(EEG.data,1)
%     meanVal = mean(EEG.data(chan,:,:), 'all'); 
%     sdVal = std(EEG.data(chan,:,:), [], 'all'); 
%     EEG.zDat(chan,:,:) = (EEG.data(chan,:,:) - meanVal) ./ sdVal; 
% end


    datMean = squeeze(std(EEG.data,[],[2,3])); 
    EEG.data(datMean==0,:,:) = []; 
    
    nbChanOrig = size(EEG.data,1); 
    nbTrialOrig = size(EEG.data,3); 

    maxDeflection = zeros(size(EEG.data,1), size(EEG.data,3)); 
    window = 80/ (1000 / EEG.srate);
    for chan = 1:size(EEG.data,1)
      
        curChan = squeeze(EEG.data(chan,:,:)); 
        deflectSlice = maxDeflection(chan, :); 
        for trial = 1:size(curChan,2)
            maxVals = movmax(curChan(:,trial), window); 
            minVals = movmin(curChan(:,trial), window);
            curDeflect = maxVals - minVals;
            deflectSlice(trial) = max(curDeflect); 


        end
        maxDeflection(chan, :) = deflectSlice; 
        
    end
    
figure
imagesc(maxDeflection)
caxis([99,100])


blinkTrials = arrayfun(@(x) sum(maxDeflection([8,14,9,21,25,22],x)>100),...
    [1:nbTrialOrig]);
%ID candidate blinks as any trials where more than 4 of the above eye 
%electrodes have large deflections
blinkTrials = find(blinkTrials>4); 
delete = []; 
for bb = 1:length(blinkTrials)
%     figure
%     plot(squeeze(EEG.data([8,14,9,21,25,22,127,126], :, blinkTrials(bb)))')
    [~, idxVals] = max(EEG.data([8,14,9,21,25,22], :, blinkTrials(bb)),[],2 );
    %they need to have their deflections within an 80ms window and have 
    %the two under eye electrodes be negative at the time of the peak to be
    %a blink
    if max(idxVals) - min(idxVals) > 80  || ... 
       sum(EEG.data([127,126], round(mean(idxVals)), blinkTrials(bb)) < 0) ~= 2
        delete = [delete, bb]; 
%         title('not a blink')
    end
    
end
blinkTrials(delete) = []; 


noiseTrials = blinkTrials; 
if ~isempty(noiseTrials)
    EEG.data(:,:,noiseTrials) = []; 
    maxDeflection(:,noiseTrials) = []; 
    behDat(noiseTrials,:) = []; 
end
    %remove channels where over 50% of trials involve a max deflection of
    %greater than 100 microvolts
    noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.30);


    if length(noiseChans) > size(EEG.data,1)/4 %if over a quarter of channels are about to be removed, then try doing rough trial removal first
        %remove trials where over 75% of channels have 100 microvolt
        %deflections
        noiseTrials = find(sum(maxDeflection>100,1)> (size(EEG.data,1)*3/4) );
        if ~isempty(noiseTrials)
            EEG.data(:,:,noiseTrials) = []; 
            maxDeflection(:,noiseTrials) = []; 
            behDat(noiseTrials,:) = [];
        end
    
        nbChanFinal = size(EEG.data,1); 
        nbTrialFinal = size(EEG.data,3);
        % then go back and do channel removal
        noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.50);
        if ~isempty(noiseChans) 
             %checking to see if there are channels with super high power
             %anywhere
             EEG.removedChans = EEG.chanlocs(noiseChans); 
             EEG.chanlocs(noiseChans) = []; 
             EEG.nbchan = EEG.nbchan - length(noiseChans); 
             EEG.data(noiseChans,:,:) = []; 
             maxDeflection(noiseChans,:) = []; 
        else
            EEG.removedChans = []; 
        end
         %remove trials where over 25% of channels have 100 microvolt
        %deflections
        noiseTrials = find(sum(maxDeflection>100,1)> (size(EEG.data,1)/4) );
        if ~isempty(noiseTrials)
            EEG.data(:,:,noiseTrials) = []; 
            maxDeflection(:,noiseTrials) = []; 
            behDat(noiseTrials,:) = []; 
        end
    
        nbChanFinal = size(EEG.data,1); 
        nbTrialFinal = size(EEG.data,3);

    else 

    
        if ~isempty(noiseChans) 
             %checking to see if there are channels with super high power
             %anywhere
             EEG.removedChans = EEG.chanlocs(noiseChans); 
             EEG.chanlocs(noiseChans) = []; 
             EEG.nbchan = EEG.nbchan - length(noiseChans); 
             EEG.data(noiseChans,:,:) = []; 
             maxDeflection(noiseChans,:) = []; 
        else
            EEG.removedChans = []; 
        end
         %remove trials where over 25% of channels have 100 microvolt
        %deflections
        noiseTrials = find(sum(maxDeflection>100,1)> (size(EEG.data,1)/4) );
        if ~isempty(noiseTrials)
            EEG.data(:,:,noiseTrials) = []; 
            maxDeflection(:,noiseTrials) = []; 
            behDat(noiseTrials,:) = []; 
        end
    
        nbChanFinal = size(EEG.data,1); 
        nbTrialFinal = size(EEG.data,3);
    end

%% rereference the data to average
   
    EEG.data = EEG.data - mean(EEG.data, 1); 

   

    EEG.noiseRemoved = true; 
    EEG.nbChanOrig = nbChanOrig; 
    EEG.nbTrialOrig = nbTrialOrig;



end