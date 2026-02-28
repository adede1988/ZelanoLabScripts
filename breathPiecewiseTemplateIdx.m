function [idx50, lm] = breathPiecewiseTemplateIdx(chanDat, varargin)
% breathPiecewiseTemplateIdx
% Returns an nBreaths x 50 matrix of indices (into the time dimension of
% chanDat.trial.rsp) that samples a standardized, piecewise "breath template":
%
%   10 pts: inhale onset  -> inhale peak                (exclude peak)
%   10 pts: inhale peak   -> first cross below onset    (exclude crossing)
%   10 pts: that crossing -> exhale trough              (exclude trough)
%   10 pts: exhale trough -> return cross above onset   (exclude crossing)
%   10 pts: pause (return crossing -> next inhale onset) (exclude next onset)
%          (pause block is all-NaN if <10 samples available)
%
% Assumptions (per user):
%   - chanDat.trial.rsp is breaths x time, spirometry velocity, inhale positive
%   - breath onset is at onsetIdx (default 1000) in each row
%   - chanDat.behDat.length is breath length (sec) from this onset to next onset
%   - data are cleaned; required fields exist
%
% Optional name/value parameters:
%   'onsetIdx'        (default 1000)
%   'nPerSeg'         (default 10)   % must be 10 for 50 total, but kept general
%   'onsetMedWinSec'  (default 0.05) % median window half-width around onset
%   'epsFrac'         (default 0.03) % epsilon = epsFrac * (peakVal - onsetLevel)
%   'minAbsEps'       (default 1e-6)
%
% Outputs:
%   idx50 : nBreaths x (5*nPerSeg) indices (NaN row if any of first 4 segments fail)
%   lm    : struct of landmark indices (useful for QC/debugging)

% -------------------- params --------------------
p = inputParser;
p.addParameter('onsetIdx',       1000,  @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('nPerSeg',        10,    @(x) isnumeric(x) && isscalar(x) && x>=1);
p.addParameter('onsetMedWinSec', 0.05,  @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('epsFrac',        0.03,  @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('minAbsEps',      1e-6,  @(x) isnumeric(x) && isscalar(x) && x>=0);
p.parse(varargin{:});
onsetIdx       = p.Results.onsetIdx;
nPerSeg        = p.Results.nPerSeg;
onsetMedWinSec = p.Results.onsetMedWinSec;
epsFrac        = p.Results.epsFrac;
minAbsEps      = p.Results.minAbsEps;

% -------------------- inputs --------------------
rsp = chanDat.trial.rsp;
fs  = chanDat.fs;

[nBreaths, nTime] = size(rsp);
if onsetIdx > nTime
    error('onsetIdx (%d) exceeds number of timepoints in rsp (%d).', onsetIdx, nTime);
end

% Extract breath length (sec) robustly
if istable(chanDat.behDat)
    try
        lenSec = chanDat.behDat.length;
    catch %need to deal with cued sniffs that don't have a length input
        
        %first derivative used to find peaks and troughs
        difVals = diff(chanDat.trial.lowRsp, [], 2);
        
        bidx = (1:nBreaths).';  % column
        
        % First peak after sample 1000: + -> - in derivative
        peaks = arrayfun(@(x) ...
            find(difVals(x,1000:end-1) > 0 & difVals(x,1001:end) < 0, 1), ...
            1:nBreaths);
        
        % First trough after the peak (still in derivative): - -> + in derivative
        troughs = arrayfun(@(x,p) ...
            find(difVals(x,p+1000:end-1) < 0 & difVals(x,p+1+1000:end) > 0, 1), ...
            bidx, peaks(:));
        
        % Convert peak/trough indices into rsp indices (offset by 1000)
        peakIdx   = 1000 + peaks(:);
        troughIdx = 1000 + peaks(:) + troughs(:);
        
        peakVals   = arrayfun(@(x,i) rsp(x,i), bidx, peakIdx);
        troughVals = arrayfun(@(x,i) rsp(x,i), bidx, troughIdx);
        
        endThresh = epsFrac * 3 + median(rsp, 'all');
     


        % First time after trough that rsp exceeds endThresh
        endIdx = arrayfun(@(x,i0) find(chanDat.trial.lowRsp(x,i0:end) > endThresh, 1), bidx, troughIdx);
        endIdx = endIdx + troughIdx;
        lenSec = (endIdx - 1000) ./ fs;
    end
elseif isstruct(chanDat.behDat)
    if numel(chanDat.behDat) > 1
        lenSec = [chanDat.behDat.length]';
    else
        lenSec = chanDat.behDat.length(:);
    end
else
    lenSec = chanDat.behDat.length(:);
end
if numel(lenSec) ~= nBreaths
    error('chanDat.behDat.length must have one value per breath (expected %d, got %d).', ...
        nBreaths, numel(lenSec));
end

% Preallocate
idx50 = nan(nBreaths, 5*nPerSeg);

lm = struct();
lm.onsetIdx      = nan(nBreaths,1); 
lm.winEnd        = nan(nBreaths,1);
lm.onsetLevelMed = nan(nBreaths,1);
lm.peakIdx       = nan(nBreaths,1);
lm.troughIdx     = nan(nBreaths,1);
lm.crossBelowIdx = nan(nBreaths,1);
lm.crossAboveIdx = nan(nBreaths,1);
lm.pauseEndIdx   = nan(nBreaths,1);
lm.epsilon       = nan(nBreaths,1);

% Median window around onset
winHalf = round(onsetMedWinSec * fs);

for b = 1:nBreaths
    x = rsp(b,:);

    % Breath end based on length (cap to available epoch)
    endIdx = onsetIdx + round(lenSec(b) * fs);
    winEnd = min(endIdx, nTime);
    lm.winEnd(b) = winEnd;

    if winEnd <= onsetIdx
        % Not enough data after onset -> invalid breath
        continue;
    end



    % Onset level (median in a short window)
    w0 = max(1, onsetIdx - winHalf);
    w1 = min(nTime, onsetIdx + winHalf);
    onsetLevel = median(x(w0:w1), 'omitnan');
    lm.onsetLevelMed(b) = onsetLevel;

    % Find inhale peak (max) within [onsetIdx, winEnd]
    w = onsetIdx:winEnd;
    [peakVal, relPeak] = max(x(w));
    peakIdx = w(relPeak);
    lm.peakIdx(b) = peakIdx;

    inhaleAmp = peakVal - onsetLevel;
    if ~(isfinite(inhaleAmp) && inhaleAmp > 0)
        % No clear positive inhale peak relative to onset level -> invalid
        continue;
    end

    epsVal = max(minAbsEps, epsFrac * inhaleAmp);
    lm.epsilon(b) = epsVal;

    % First crossing BELOW onsetLevel-eps after peak (within breath window)
    cbSearchStart = peakIdx;
    cbSearchEnd   = winEnd;
    cb = find(x(cbSearchStart:cbSearchEnd) <= (onsetLevel - epsVal), 1, 'first');
    if isempty(cb)
        continue;
    end
    crossBelowIdx = cbSearchStart + cb - 1;
    lm.crossBelowIdx(b) = crossBelowIdx;

    % Exhale trough (min) AFTER the below-crossing (keeps correct ordering)
    trSearchStart = crossBelowIdx;
    trSearchEnd   = winEnd;
    [~, relTr] = min(x(trSearchStart:trSearchEnd));
    troughIdx = trSearchStart + relTr - 1;
    lm.troughIdx(b) = troughIdx;

    % Return crossing ABOVE onsetLevel+eps after trough (end of exhale fall)
    caSearchStart = troughIdx;
    caSearchEnd   = winEnd;

    ca = find(x(caSearchStart:caSearchEnd) >= (median(x) - epsVal), 1, 'first');
    if isempty(ca)
        continue;
    end
    crossAboveIdx = caSearchStart + ca - 1;
    lm.crossAboveIdx(b) = crossAboveIdx;

    % Pause end: next inhale onset index based on breath length.
    % If next onset lies beyond epoch, treat pause as "to end of epoch".
    pauseEnd = onsetIdx + round(lenSec(b) * fs);
    pauseEnd = min(pauseEnd, nTime + 1);  % allow "to end" by using nTime+1 as exclusive bound
    lm.pauseEndIdx(b) = pauseEnd;

    % Build the 5 segments (include first point, drop last point)
    seg1 = sampleIdx(onsetIdx,      peakIdx,       nPerSeg);
    seg2 = sampleIdx(peakIdx,       crossBelowIdx, nPerSeg);
    seg3 = sampleIdx(crossBelowIdx, troughIdx,     nPerSeg);
    seg4 = sampleIdx(troughIdx,     crossAboveIdx, nPerSeg);

    % First 4 segments are mandatory (all-or-nothing breath validity)
    if isempty(seg1) || isempty(seg2) || isempty(seg3) || isempty(seg4)
        continue;
    end

    % Pause segment is optional: all-NaN if insufficient samples
    seg5 = sampleIdx(crossAboveIdx, pauseEnd, nPerSeg);
    if isempty(seg5)
        seg5 = nan(1, nPerSeg);
    end

    idx50(b,:) = [seg1, seg2, seg3, seg4, seg5];
end



%% go through again and correct NaNs 

binnedVals = [-200:5:200];
binnedCounts = arrayfun(@(x) sum(lm.onsetLevelMed<x), binnedVals);
binnedCounts = diff(binnedCounts); 

[~, onsetMaxIdx] = max(binnedCounts);  

lowestValid = find(flip(binnedCounts(1:onsetMaxIdx))/...
                                        sum(binnedCounts)<.05,1);
lowestValid = onsetMaxIdx - lowestValid; 

highestValid = find(binnedCounts(onsetMaxIdx:end)/...
                                        sum(binnedCounts)<.05,1); 
highestValid = highestValid + onsetMaxIdx; 

lm.onsetLevelMed(lm.onsetLevelMed<binnedVals(lowestValid)) = ...
                                                  binnedVals(lowestValid); 
lm.onsetLevelMed(lm.onsetLevelMed>binnedVals(highestValid)) = ...
                                                  binnedVals(highestValid);

for b = 1:nBreaths
    x = rsp(b,:);

    % Breath end based on length (cap to available epoch)
    endIdx = onsetIdx + round(lenSec(b) * fs);
    winEnd = min(endIdx, nTime);
    lm.winEnd(b) = winEnd;

    if winEnd <= onsetIdx
        % Not enough data after onset -> invalid breath
        continue;
    end



    % Onset level (median in a short window)
    w0 = max(1, onsetIdx - winHalf);
    w1 = min(nTime, onsetIdx + winHalf);
    onsetLevel = median(x(w0:w1), 'omitnan');
    if onsetLevel ~= lm.onsetLevelMed(b)
        onsetLevel = lm.onsetLevelMed(b);
        onsetCurIdx = onsetIdx - 150 + ...
                        find(x(onsetIdx-150:onsetIdx+100)>=onsetLevel, 1);
        if isempty(onsetCurIdx)
            idx50(b,:) = nan(1, nPerSeg*5);
            continue;
        end
        lm.onsetIdx(b) = onsetCurIdx; 
    else
        lm.onsetIdx(b) = onsetIdx; 
        onsetCurIdx = onsetIdx; 
    end
        
    % Find inhale peak (max) within [onsetIdx, winEnd]
    w = onsetCurIdx:winEnd;
    [peakVal, relPeak] = max(x(w));
    peakIdx = w(relPeak);
    lm.peakIdx(b) = peakIdx;

    inhaleAmp = peakVal - onsetLevel;
    if ~(isfinite(inhaleAmp) && inhaleAmp > 0)
        % No clear positive inhale peak relative to onset level -> invalid
        continue;
    end

    epsVal = max(minAbsEps, epsFrac * inhaleAmp);
    lm.epsilon(b) = epsVal;

    % First crossing BELOW onsetLevel-eps after peak (within breath window)
    cbSearchStart = peakIdx;
    cbSearchEnd   = winEnd;
    cb = find(x(cbSearchStart:cbSearchEnd) <= (onsetLevel - epsVal), 1, 'first');
    if isempty(cb)
        continue;
    end
    crossBelowIdx = cbSearchStart + cb - 1;
    lm.crossBelowIdx(b) = crossBelowIdx;

    % Exhale trough (min) AFTER the below-crossing (keeps correct ordering)
    trSearchStart = crossBelowIdx;
    trSearchEnd   = winEnd;
    [~, relTr] = min(x(trSearchStart:trSearchEnd));
    troughIdx = trSearchStart + relTr - 1;
    lm.troughIdx(b) = troughIdx;

    % Return crossing ABOVE onsetLevel+eps after trough (end of exhale fall)
    caSearchStart = troughIdx;
    caSearchEnd   = winEnd;
    ca = find(x(caSearchStart:caSearchEnd) >= (onsetLevel - epsVal), 1, 'first');
    if isempty(ca)
        continue;
    end
    crossAboveIdx = caSearchStart + ca - 1;
    lm.crossAboveIdx(b) = crossAboveIdx;

    % Pause end: next inhale onset index based on breath length.
    % If next onset lies beyond epoch, treat pause as "to end of epoch".
    pauseEnd = onsetIdx + round(lenSec(b) * fs);
    pauseEnd = min(pauseEnd, nTime + 1);  % allow "to end" by using nTime+1 as exclusive bound
    lm.pauseEndIdx(b) = pauseEnd;

    % Build the 5 segments (include first point, drop last point)
    seg1 = sampleIdx(onsetCurIdx,      peakIdx,       nPerSeg);
    seg2 = sampleIdx(peakIdx,       crossBelowIdx, nPerSeg);
    seg3 = sampleIdx(crossBelowIdx, troughIdx,     nPerSeg);
    seg4 = sampleIdx(troughIdx,     crossAboveIdx, nPerSeg);

    % First 4 segments are mandatory (all-or-nothing breath validity)
    if isempty(seg1) || isempty(seg2) || isempty(seg3) || isempty(seg4)
        continue;
    end

    % Pause segment is optional: all-NaN if insufficient samples
    seg5 = sampleIdx(crossAboveIdx, pauseEnd, nPerSeg);
    if isempty(seg5)
        seg5 = nan(1, nPerSeg);
    end

    idx50(b,:) = [seg1, seg2, seg3, seg4, seg5];
end



end % function


% -------------------------------------------------------------------------
function seg = sampleIdx(a, b, n)
% Return 1 x n integer indices sampling [a, b) linearly
% include a, exclude b (drop last point)
% Requires at least n samples available in [a, b) => (b - a) >= n
seg = [];
if ~isfinite(a) || ~isfinite(b) || b <= a
    return;
end
a = round(a); b = round(b);
if (b - a) < n
    return;
end
% Monotone, integer, unique by construction under (b-a)>=n
seg = a + floor((0:n-1) * (b - a) / n);
seg = double(seg(:))'; % row
end
