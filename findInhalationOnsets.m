function [onsets, diagnostics] = findInhalationOnsets(rspSig, fs, peaks, troughs)
% FINDINHALATIONONSETS
%   Detect inhalation onset times from a respiration velocity signal.
%
%   INPUTS
%     rspSig   : respiration velocity time series (column vector)
%     fs       : sampling rate (Hz)
%     peaks    : indices of inhalation peaks (one per breath)
%     troughs  : indices of preceding troughs (one per breath)
%
%   OUTPUTS
%     onsets      : index of inferred inhalation onset for each breath
%     diagnostics : N x 3 matrix with per-breath diagnostics:
%                    col 1 = rspSig value at final chosen onset
%                    col 2 = index of final chosen onset
%                    col 3 = final "choice" value (biased difVal score)
%
%   The algorithm has two passes:
%     1) For each breath, estimate an initial onset between the trough and
%        peak using local slope and MSE fits.
%     2) Refine onset by scanning candidate points before the peak and
%        picking the location that maximizes a slope-difference metric
%        (left vs right) *biased* toward mid-range amplitudes.
%
%   NOTE: This is a direct packaging of the provided script into a
%   function; no changes to the underlying logic have been made.

% downsample data for time
smthRsp = rspSig; 
rspSig = rspSig(1:10:end); 
fs = round(fs/10); 
peaks = round(peaks/10); 
troughs = round(troughs/10); 


% -------------------------------------------------------------------------
% FIRST PASS: Initial onset estimate between trough and peak
% -------------------------------------------------------------------------
onsets        = zeros(size(peaks)); 
onsetHeights  = onsets; 

for ii = 1:length(peaks)

    % For each breath, search between the peak and the preceding trough
    % to find the initial inhale onset.
    ti = troughs(ii); 
    pi = peaks(ii); 

    % Preallocate slope and MSE arrays for candidate onset locations
    slopes = zeros(pi - ti - round(fs/40), 1);
    MSE    = slopes; 

    % For each candidate jj between trough and (peak - small window),
    % compute slope of line between smthRsp(jj) and smthRsp(pi),
    % and mean squared error of that line fit.
    for jj = ti : pi - round(fs/40)
        slopes(jj - ti + 1) = (smthRsp(pi) - smthRsp(jj)) / (pi - jj); 
        
        fittedVals = linspace(smthRsp(jj), smthRsp(pi), (pi - jj + 1)); 
        MSE(jj - ti + 1) = sum((fittedVals - smthRsp(jj:pi)).^2) / (pi - jj); 
    end

    % Find an initial inflection as the max slope
    [~, inflection] = max(slopes(1:end - round(fs/40))); 

    % Use MSE distribution to refine inflection index
    meanValLocal = mean(MSE(inflection:end));
    sdValLocal   = std(MSE(inflection:end));

    inflection = length(MSE) - find(flip(MSE) < meanValLocal + sdValLocal, 1, 'last');
    onsets(ii)       = inflection + ti; 
    onsetHeights(ii) = rspSig(inflection + ti); 
end

% -------------------------------------------------------------------------
% SECOND PASS SETUP: Use distribution of first-pass onsets as bias
% -------------------------------------------------------------------------
meanVal = mean(onsetHeights); 
sdVal   = std(onsetHeights); %#ok<NASGU>  % kept for parity with original code

% Mean trough and peak values (used for amplitude biasing)
meanT = mean(rspSig(troughs)); 
meanP = mean(rspSig(peaks)); 

% Diagnostics: [onsetAmplitude, onsetIndex, choiceValue]
diagnostics = zeros(length(onsets), 3); 

% -------------------------------------------------------------------------
% SECOND PASS: Refine onset using left vs right slopes + amplitude bias
% -------------------------------------------------------------------------
parfor ii = 1:length(peaks) 
    % Local diagnostics slice for this breath
    slice = diagnostics(ii,:); 
    
    ti = troughs(ii); 
    pi = peaks(ii); 

    % Window for right-side slope estimate
    win = min(round(fs/3), round((pi - ti) / 5));

    % Candidates: search from (peak - up to 2 seconds) back to trough
    searchStart = max(pi - fs*2, ti); 
    candidates  = searchStart:pi; 

    % Prepare a padded test signal to avoid out-of-range when indexing
    testSig = smthRsp; 
    testSig(pi : pi + fs) = testSig(pi); 
    testSig(ti - fs*2 : ti) = testSig(ti); 

    % Steps at which to evaluate left slopes
    steps = round(linspace(round(fs/8), win, 5)); 

    difVal = zeros(length(candidates), 1); 

    % If no candidates (degenerate case), fall back to mid-point
    if isempty(candidates)
        onsets(ii)        = round((pi + ti) / 2);
        slice(2)          = round((pi + ti) / 2); 
        diagnostics(ii,:) = slice; 
        continue
    end

    % ---------------------------------------------------------------------
    % For each candidate onset:
    %   - Compute left-side slopes over multiple window lengths
    %   - Compute right-side slope over fixed window
    %   - Store mean(rightSlope - leftSlopes) in difVal
    % ---------------------------------------------------------------------
    for jj = 1:length(candidates)
        idx = candidates(jj); 

        % ---- Left side slopes (multiple window scales) ----
        yLeft = arrayfun(@(x) testSig(idx - x : idx), steps, ...
                         'UniformOutput', false);
        xLeft = cellfun(@(yy) 1:length(yy), yLeft, 'UniformOutput', false);
        rLeft = cellfun(@(xx,yy) corr(xx(:), yy(:)), xLeft, yLeft, ...
                        'UniformOutput', false);
        leftSlopes = cellfun(@(xx,yy,zz) zz * std(yy)/std(xx), ...
                             xLeft, yLeft, rLeft);

        % ---- Right side slope (single window) ----
        ridx = idx + win; 
        yRight = testSig(idx:ridx);
        xRight = 1:length(yRight);
        rRight = corr(xRight(:), yRight(:));
        rightSlope = rRight * std(yRight)/std(xRight);

        % Store slope difference metric
        difVal(jj) = mean(rightSlope - leftSlopes); 
    end

    % ---------------------------------------------------------------------
    % Amplitude biasing:
    %   - Prefer candidates near meanVal
    %   - Suppress candidates near mean peak (meanP) or trough (meanT)
    % ---------------------------------------------------------------------
    amp = rspSig(candidates);   % respiration velocity at candidate points
    
    % Preallocate amplitude weight
    wAmp = zeros(size(amp));
    
    % Above the mean (between meanVal and meanP)
    idxUp = amp >= meanVal;
    if any(idxUp)
        fracUp = (amp(idxUp) - meanVal) ./ (meanP - meanVal);   % 0 at mean, 1 at peak
        wAmp(idxUp) = 1 - fracUp.^2;                            % 1 at mean, 0 at peak
    end
    
    % Below the mean (between meanT and meanVal)
    idxDown = amp < meanVal;
    if any(idxDown)
        fracDown = (meanVal - amp(idxDown)) ./ (meanVal - meanT); % 0 at mean, 1 at trough
        wAmp(idxDown) = 1 - fracDown.^2;                          % 1 at mean, 0 at trough
    end
    
    % Clamp in case of slight numerical overshoot
    wAmp = max(wAmp, 0);
    
    % Combine bias with original detection strength
    test = difVal(:) .* wAmp(:);

    % ---------------------------------------------------------------------
    % Final choice for this breath: max biased score
    % ---------------------------------------------------------------------
    [choiceVal, idxChoice] = max(test); 
    
    onsets(ii) = candidates(idxChoice); 

    slice(1) = rspSig(candidates(idxChoice));  % amplitude at chosen onset
    slice(2) = onsets(ii);                     % onset index
    slice(3) = choiceVal;                      % biased score
    diagnostics(ii,:) = slice; 
end


onsets = onsets*10; 

end
