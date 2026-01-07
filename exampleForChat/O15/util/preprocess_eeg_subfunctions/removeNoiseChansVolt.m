function [badTS, badChans] = removeNoiseChansVolt(EEG, fs)

%input: 
%       EEG             channels X time matrix of EEG data
%       labels          1 X channels cell array of channel labels
%       fs              sampling rate in Hz

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
%Fall 2025

%split data into 2 second trials

[c, L] = size(EEG); 

snipL = fs*2; 
starts = [1:snipL:L]; 
tmpEEG = zeros([c,snipL,length(starts)]); 

for ii = 1:length(starts)
    snipStart = starts(ii); 
    snipEnd = min(L, snipStart+snipL-1);
    tmpL = snipEnd - snipStart + 1; 
    tmpEEG(:,1:tmpL,ii) = EEG(:,snipStart:snipEnd); 
    
end


maxDeflection = zeros(size(tmpEEG,1), size(tmpEEG,3)); 
badRecord = zeros(size(maxDeflection)); 
window = round(80/ (1000 / fs));
for chan = 1:size(tmpEEG,1)
    for trial = 1:size(tmpEEG,3)
        for ii = 1:size(tmpEEG,2) - window
            cur = max(tmpEEG(chan,ii:ii+window, trial)) - min(tmpEEG(chan,ii:ii+window, trial));
            if cur > maxDeflection(chan, trial)
                maxDeflection(chan,trial) = cur; 
            end
        end
    end
end
    

%remove channels where over 50% of trials involve a max deflection of
%greater than 100 microvolts
noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.50);


if length(noiseChans) > size(tmpEEG,1)/4 %if over a quarter of channels are about to be removed, then try doing rough trial removal first
    %remove trials where over 75% of channels have 100 microvolt
    %deflections
    noiseTrials = find(sum(maxDeflection>100,1)> (size(tmpEEG,1)*3/4) );
    if ~isempty(noiseTrials)
        badRecord(:,noiseTrials) = 1; 
        maxDeflection(:,noiseTrials) = 0; 
    end

   
    % then go back and do channel removal
    noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.50);
    if ~isempty(noiseChans) 
        badRecord(noiseChans,:) = 1; 
        maxDeflection(noiseChans,:) = 0; 
    end
     %remove trials where over 25% of channels have 100 microvolt
    %deflections
    noiseTrials = find(sum(maxDeflection>100,1)> (size(tmpEEG,1)/4) );
    if ~isempty(noiseTrials)
       badRecord(:,noiseTrials) = 1; 
       maxDeflection(:,noiseTrials) = 0; 
    end


else 


    if ~isempty(noiseChans) 
        badRecord(noiseChans,:) = 1; 
        maxDeflection(noiseChans,:) = 0; 
    end
     %remove trials where over 25% of channels have 100 microvolt
    %deflections
    noiseTrials = find(sum(maxDeflection>100,1)> (size(tmpEEG,1)/4) );
    if ~isempty(noiseTrials)
       badRecord(:,noiseTrials) = 1; 
       maxDeflection(:,noiseTrials) = 0; 
    end

  
end

badChans = find(sum(badRecord,2)==length(starts));
badTrials = sum(badRecord,1)==c;
badTrials(end) = false; 
badTS = zeros(L,1); 
for ii = 1:length(starts)
    if badTrials(ii)
        badTS(starts(ii):starts(ii)+snipL-1) = 1; 
    end
end
   

end
