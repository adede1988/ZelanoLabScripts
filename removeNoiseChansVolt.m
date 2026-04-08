function [badTS, badChans, newEEG, interpChan] = removeNoiseChansVolt(EEG, fs, skipChans, chanLocs)

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

% 
% maxDeflection = zeros(size(tmpEEG,1), size(tmpEEG,3)); 
% badRecord = zeros(size(maxDeflection)); 
% window = round(5/ (1000 / fs));
% for chan = 1:size(tmpEEG,1)
%     for trial = 1:size(tmpEEG,3)
%         for ii = 1:round(window/4):(size(tmpEEG,2) - window)
%             cur = max(tmpEEG(chan,ii:ii+window, trial)) - min(tmpEEG(chan,ii:ii+window, trial));
%             if cur > maxDeflection(chan, trial)
%                 maxDeflection(chan,trial) = cur; 
%             end
%         end
%     end
% end




%%% chat speed up: 

[nChan, nTime, nTrial] = size(tmpEEG);

window = round(10 / (1000 / fs));   % same as your code
step   = round(window / 4);

startIdx = 1:step:(nTime - window);

% reshape to 2D: each row = one chan/trial
x = permute(tmpEEG, [1 3 2]);       % [chan x trial x time]
x = reshape(x, [], nTime);          % [(chan*trial) x time]

% for each start point i, compute max/min over i:(i+window)
winRange = movmax(x, [0 window], 2) - movmin(x, [0 window], 2);

% keep only the sampled window starts, then take the largest deflection
maxDeflection = max(winRange(:, startIdx), [], 2);

% reshape back to [chan x trial]
maxDeflection = reshape(maxDeflection, nChan, nTrial);

badRecord = zeros(size(maxDeflection));












  
%%%%%%%%%%START NEW
noiseTrials = find(sum(maxDeflection>100,1)> (size(tmpEEG,1)*3/4) ); %first remove 75% bad trials entirely, these are unusable
if ~isempty(noiseTrials)
    badRecord(:,noiseTrials) = 1; 
    maxDeflection(:,noiseTrials) = 0; 
end
%do blockwise noise channel detection
blockSize = 5; 
badChans = zeros(size(maxDeflection)); 
for block = 1:size(maxDeflection,2)-blockSize
    curNoiseChans = find(sum(maxDeflection(:,block:block+blockSize)>50,2) > (blockSize/2) );
    curNoiseChans2 = find(sum(maxDeflection(:,block:block+blockSize)<5,2) > (blockSize/2) );
    idx = ismember(curNoiseChans, skipChans);
    if ~isempty(idx)
        curNoiseChans(idx) = []; 
    end
    badChans(curNoiseChans,block:block+blockSize) = 1; 
    badChans(curNoiseChans2,block:block+blockSize) = 1; 
end
outBad = find(sum(badChans,2) > size(tmpEEG,3)*.5);
interpCount = zeros(size(tmpEEG, 3),1); 
for tt = 1:size(tmpEEG,3)
    if(sum(badChans(:,tt)) > 0)
        if sum(badChans(:,tt)) < (nChan/2)
        tmpEEG(:,:,tt) = interpolate_perrinX(tmpEEG(:,:,tt), ...
                                                 chanLocs.X, ...
                                                 chanLocs.Y, ...
                                                 chanLocs.Z, ...
                                                 find(badChans(:,tt)==1));
        interpCount(tt) = sum(badChans(:,tt));
            
        else
            badRecord(:,badChans(:,tt)==1) = 1; 
        end
    end
end

newEEG = zeros(size(EEG)); 
interpChan = zeros(size(newEEG,2),1);
for tt = 1:size(tmpEEG,3)
    try
    newEEG(:, starts(tt):min(L,starts(tt)+snipL-1) ) = tmpEEG(:,:,tt); 
    interpChan(starts(tt):min(L,starts(tt)+snipL-1)) = interpCount(tt); 
    catch
    end
end

figure; 
hold on 
for ii = 1:32

    plot(EEG(ii,:) + ii*50, 'color', 'green')
end

for ii = 1:32

    plot(newEEG(ii,:) + ii*50, 'color', 'k')
end



%% recalculate the max deflections: 

snipL = fs*2; 
starts = [1:snipL:L]; 
tmpEEG = zeros([c,snipL,length(starts)]); 

for ii = 1:length(starts)
    snipStart = starts(ii); 
    snipEnd = min(L, snipStart+snipL-1);
    tmpL = snipEnd - snipStart + 1; 
    tmpEEG(:,1:tmpL,ii) = newEEG(:,snipStart:snipEnd); 
    
end

% 
% maxDeflection = zeros(size(tmpEEG,1), size(tmpEEG,3)); 
% badRecord = zeros(size(maxDeflection)); 
% window = round(5/ (1000 / fs));
% for chan = 1:size(tmpEEG,1)
%     for trial = 1:size(tmpEEG,3)
%         for ii = 1:round(window/4):(size(tmpEEG,2) - window)
%             cur = max(tmpEEG(chan,ii:ii+window, trial)) - min(tmpEEG(chan,ii:ii+window, trial));
%             if cur > maxDeflection(chan, trial)
%                 maxDeflection(chan,trial) = cur; 
%             end
%         end
%     end
% end




%%% chat speed up: 

[nChan, nTime, nTrial] = size(tmpEEG);

window = round(80 / (1000 / fs));   % same as your code
step   = round(window / 4);

startIdx = 1:step:(nTime - window);

% reshape to 2D: each row = one chan/trial
x = permute(tmpEEG, [1 3 2]);       % [chan x trial x time]
x = reshape(x, [], nTime);          % [(chan*trial) x time]

% for each start point i, compute max/min over i:(i+window)
winRange = movmax(x, [0 window], 2) - movmin(x, [0 window], 2);

% keep only the sampled window starts, then take the largest deflection
maxDeflection = max(winRange(:, startIdx), [], 2);

% reshape back to [chan x trial]
maxDeflection = reshape(maxDeflection, nChan, nTrial);







%remove channels where over 50% of trials involve a max deflection of
%greater than 100 microvolts
noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.50);
idx = ismember(noiseChans, skipChans);
    if ~isempty(idx)
        noiseChans(idx) = []; 
    end

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

badChans = find(sum(badRecord,2)>=length(starts)*.6);
badTrials = sum(badRecord,1)==c;
badTrials(end) = false; 
badTS = zeros(L,1); 
for ii = 1:length(starts)
    if badTrials(ii)
        badTS(starts(ii):starts(ii)+snipL-1) = 1; 
    end
end

badChans = badChans(:); 
outBad = outBad(:); 
badChans = [badChans; outBad];





















%%%%%%%%%%%%END NEW








% 
% 
% %remove channels where over 50% of trials involve a max deflection of
% %greater than 100 microvolts
% noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.50);
% 
% 
% if length(noiseChans) > size(tmpEEG,1)/4 %if over a quarter of channels are about to be removed, then try doing rough trial removal first
%     %remove trials where over 75% of channels have 100 microvolt
%     %deflections
%     noiseTrials = find(sum(maxDeflection>100,1)> (size(tmpEEG,1)*3/4) );
%     if ~isempty(noiseTrials)
%         badRecord(:,noiseTrials) = 1; 
%         maxDeflection(:,noiseTrials) = 0; 
%     end
% 
% 
%     % then go back and do channel removal
%     noiseChans = find(sum(maxDeflection>100,2) ./ size(maxDeflection,2)>.50);
%     if ~isempty(noiseChans) 
%         badRecord(noiseChans,:) = 1; 
%         maxDeflection(noiseChans,:) = 0; 
%     end
%      %remove trials where over 25% of channels have 100 microvolt
%     %deflections
%     noiseTrials = find(sum(maxDeflection>100,1)> (size(tmpEEG,1)/4) );
%     if ~isempty(noiseTrials)
%        badRecord(:,noiseTrials) = 1; 
%        maxDeflection(:,noiseTrials) = 0; 
%     end
% 
% 
% else 
% 
% 
%     if ~isempty(noiseChans) 
%         badRecord(noiseChans,:) = 1; 
%         maxDeflection(noiseChans,:) = 0; 
%     end
%      %remove trials where over 25% of channels have 100 microvolt
%     %deflections
%     noiseTrials = find(sum(maxDeflection>100,1)> (size(tmpEEG,1)/4) );
%     if ~isempty(noiseTrials)
%        badRecord(:,noiseTrials) = 1; 
%        maxDeflection(:,noiseTrials) = 0; 
%     end
% 
% 
% end
% 
% badChans = find(sum(badRecord,2)==length(starts));
% badTrials = sum(badRecord,1)==c;
% badTrials(end) = false; 
% badTS = zeros(L,1); 
% for ii = 1:length(starts)
%     if badTrials(ii)
%         badTS(starts(ii):starts(ii)+snipL-1) = 1; 
%     end
% end
   

end
