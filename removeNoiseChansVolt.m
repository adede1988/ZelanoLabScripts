function [badTS, badChans, newEEG, interpChan] = removeNoiseChansVolt(EEG, fs, skipChans, chanLocs)

%input:
%       EEG         channels X time matrix of EEG data (microvolts)
%       fs          sampling rate in Hz
%       skipChans   vector of channel indices to protect from high-amplitude
%                   flagging (e.g. EOG, which is legitimately high-amplitude
%                   and should not be removed for exceeding the voltage
%                   thresholds). These channels ARE still eligible to be
%                   flagged and removed if flat/dead.
%       chanLocs    struct with fields .X .Y .Z giving channel coordinates,
%                   used for Perrin spherical-spline interpolation

%output:
%       badTS       L X 1 binary time mask. A timepoint is 1 only when its
%                   parent 2 s trial had ALL channels rejected, i.e. trials
%                   that tripped the trial-level thresholds. Use as a
%                   "drop these segments entirely" mask. The final trial is
%                   forced to 0.
%       badChans    vector of unique channel indices to drop entirely.
%                   Combines (a) channels flagged in >=60% of trials in the
%                   second pass and (b) channels flagged in >50% of trials by
%                   the first-pass blockwise detector. These channels are NOT
%                   removed from newEEG; removal is left to the caller.
%       newEEG      channels X time cleaned data, same shape as EEG, with
%                   transient bad channels spatially interpolated per 2 s
%                   trial (only where fewer than half the channels were bad
%                   in that trial).
%       interpChan  L X 1 count, per timepoint, of how many channels were
%                   interpolated in the covering trial. A repair-intensity
%                   trace: high stretches were heavily reconstructed.

%PROCESSING OVERVIEW
%The function epochs the data into non-overlapping 2 s trials. The final
%short epoch is handled at its true length: internally it is edge-padded
%(its last real sample is replicated) only so the sliding peak-to-peak
%stats reflect its real data, and on write-back only its real samples are
%used. Two passes are then run.
%
%PASS 1 (detect + repair, ~10 ms sliding window):
%  - Compute peak-to-peak amplitude in a sliding ~10 ms window; reduce to the
%    single largest deflection per channel/trial.
%  - Reject the worst trials first: trials where >75% of channels exceed
%    100 uV are marked bad and removed from the channel statistics.
%  - Blockwise channel detection over a sliding 6-trial block: in each block,
%    flag a channel as noisy if it exceeds 50 uV in >=3 trials, or as flat if
%    it falls below 5 uV in >=3 trials. skipChans are exempt from the noisy
%    (high-amplitude) test so that EOG is preserved, but are NOT exempt from
%    the flat test. This yields a time-resolved bad-channel map (a channel
%    may be bad in some epochs and fine in others).
%  - Channels flagged in >50% of trials are held aside (outBad).
%  - For each trial with flagged channels: if fewer than half the channels
%    are bad, interpolate them with interpolate_perrinX and log the count;
%    otherwise mark the whole trial bad.
%  - Reassemble the interpolated trials into newEEG and build interpChan.
%  - Plot first 32 channels (green = original, black = cleaned) for QC.
%
%PASS 2 (decide what is irreparable, ~80 ms sliding window):
%  - Re-epoch newEEG and recompute max deflection with an ~80 ms window
%    (sensitive to slower, larger swings rather than fast spikes).
%  - Flag channels exceeding 100 uV in >50% of trials (skipChans exempt, to
%    preserve EOG -- including the heavy-noise recompute branch).
%  - If that would remove >25% of channels, do a rough trial purge
%    (>75% channels bad), recompute channels, then a moderate trial purge
%    (>25% channels bad). Otherwise remove the flagged channels, then purge
%    trials with >25% channels bad.
%  - Collapse the accumulated bad-record into badTS, badChans, and outBad.

%window sizes (10 ms / 80 ms) and the inline comments below are the live
%values; the commented-out loop versions are legacy and may differ.

%written by Adam Dede (adam.osman.dede@gmail.com)
%Fall 2025

%split data into 2 second trials

[c, L] = size(EEG); 

snipL  = fs*2; 
starts = 1:snipL:L; 
nTrial = length(starts); 
tmpEEG   = zeros([c, snipL, nTrial]); 
epochLen = zeros(nTrial, 1);   % true (unpadded) length of each epoch

for ii = 1:nTrial
    snipStart = starts(ii); 
    snipEnd   = min(L, snipStart+snipL-1);
    tmpL      = snipEnd - snipStart + 1; 
    epochLen(ii) = tmpL; 
    tmpEEG(:,1:tmpL,ii) = EEG(:,snipStart:snipEnd); 
    if tmpL < snipL
        % Final short epoch: pad the remainder with the last real sample
        % (edge replication) instead of zeros. The padded region then
        % contributes ~0 to the sliding peak-to-peak, so the deflection
        % stats reflect only the real data. The epoch is written back at
        % its true length below.
        tmpEEG(:, tmpL+1:snipL, ii) = repmat(EEG(:,snipEnd), 1, snipL - tmpL); 
    end
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
            badRecord(badChans(:,tt)==1,tt) = 1; 
        end
    end
end

newEEG     = zeros(size(EEG)); 
interpChan = zeros(size(newEEG,2),1);
for tt = 1:nTrial
    vLen = epochLen(tt); 
    cols = starts(tt):starts(tt)+vLen-1;        % true-length range, never exceeds L
    newEEG(:, cols)  = tmpEEG(:, 1:vLen, tt);   % write only the real samples
    interpChan(cols) = interpCount(tt); 
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

snipL  = fs*2; 
starts = 1:snipL:L; 
nTrial = length(starts); 
tmpEEG = zeros([c, snipL, nTrial]); 

for ii = 1:nTrial
    snipStart = starts(ii); 
    snipEnd   = min(L, snipStart+snipL-1);
    tmpL      = snipEnd - snipStart + 1; 
    tmpEEG(:,1:tmpL,ii) = newEEG(:,snipStart:snipEnd); 
    if tmpL < snipL
        % same edge-padding as Pass 1 for the final short epoch
        tmpEEG(:, tmpL+1:snipL, ii) = repmat(newEEG(:,snipEnd), 1, snipL - tmpL); 
    end
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
    idx = ismember(noiseChans, skipChans);   % protect EOG from amplitude rejection here too
    if ~isempty(idx)
        noiseChans(idx) = []; 
    end
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
        badTS(starts(ii):min(L, starts(ii)+snipL-1)) = 1; 
    end
end

badChans = badChans(:); 
outBad = outBad(:); 
badChans = [badChans; outBad];

badChans = unique(badChans); 



















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
