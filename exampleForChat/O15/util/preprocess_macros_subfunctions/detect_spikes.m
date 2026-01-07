function [flags, prominence] = detect_spikes(data, threshold, ...
                    winSize, markWholeWindow, gammaDat)
% DETECT_SPIKES detect sharp spike-like deflections in 3D data
% 
%   Inputs:
%     data            - C x T x N  (Channels x Time x Trials)
%     threshold       - scalar or Cx1 vector giving the threshold for (max-min)
%     winSize         - odd integer window length (default = 5)
%     markWholeWindow - logical, if true mark all samples in the window when
%                       detection happens; if false mark only the window center (default = false)
%     requireFullWindow - logical, if true ignore edges where full centered window isn't available (default = true)
%
%   Outputs:
%     flags    - binary matrix C x T x N, 1 = spike detected (at that time sample)
%     rangeVals - the computed moving (max - min) values, same size as data
%
% Example:
%   [f,r] = detect_spikes(myData, 50, 5);
%

% Defaults
if nargin < 3 || isempty(winSize),       winSize = 5; end
if nargin < 4 || isempty(markWholeWindow), markWholeWindow = false; end
if nargin < 5 || isempty(gammaDat),      gammaDat = data; end
if mod(winSize,2)==0, error('winSize must be odd.'); end

% Accept 2D (C x T) or 3D (C x T x N)
origIs2D = ismatrix(data);
if origIs2D
    data    = reshape(data,    size(data,1), size(data,2), 1);
    gammaDat= reshape(gammaDat,size(gammaDat,1), size(gammaDat,2), 1);
end

[C, T, N] = size(data);

% accept scalar threshold or per-channel
if isscalar(threshold)
    thr = threshold;
else
    threshold = reshape(threshold, [], 1);
    if numel(threshold) ~= C
        error('If threshold is a vector it must have length equal to number of channels (C).')
    end
    thr = permute(repmat(threshold,1,T,N), [1 2 3]); % C x T x N
end

% ---- compute moving max - min with centered window ----
half = floor((winSize-1)/2);
left = half;
right = winSize - 1 - left; % equals half for odd winSize

% Use movmax / movmin with centered window [left right] along time (dim=2)
rangeVals = movmax(data, [left right], 2) - movmin(data, [left right], 2);

% ---- detection ----
detections = rangeVals > thr;
prominence = zeros(size(detections));  
%check if there is a sustained gamma oscillation
%if yes, keep

chGlobalMAD = zeros(size(gammaDat,1),1);
for ch = 1:size(gammaDat,1)
    x = reshape(gammaDat(ch,:,:), 1, []);
    chGlobalMAD(ch) = mad(x, 0);  % MATLAB's MAD; already ~1.4826*median|...|
end
flankHalf = 50;       % your ~100-sample total flank width
madFloorFrac = 0.3;   % 30% of global MAD as a safety floor

parfor ch = 1:C
    ch
    limVals = prctile(real(gammaDat(ch,:,:)), [10, 90], [2, 3]);
    lower = limVals(1); upper = limVals(2);

    for tr = 1:N

        sig = squeeze(gammaDat(ch,:,tr));  % 1 x T
        for k = 1:T
            if detections(ch,k,tr) && ...
               (gammaDat(ch,k,tr) < lower || gammaDat(ch,k,tr) > upper)
                % --------- POSITIVE spike ---------
                if sig(k) > upper
                    % first index to the RIGHT that falls back below 'upper'
                    right = find(sig(k:end) < upper, 1, 'first');
                    % make absolute
                    if ~isempty(right), right = k + right - 1; end  
                    % last index to the LEFT that is below 'upper'
                    left  = find(sig(1:k) < upper, 1, 'last');
                    % cenMax = NaN; 
                    % Lmax = NaN;
                    % Rmax = NaN; 
                    if ~isempty(right) && ~isempty(left)
                        % LL = max(1, left - flankHalf);
                        % Lmax = max(sig(LL:left));
                        % RR = min(T, right + flankHalf);
                        % Rmax = max(sig(right:RR));
                        % cenMax = max(sig(left:right));
                        % prominence(ch,k,tr) = cenMax / max([Lmax, Rmax, eps]);

                    elseif ~isempty(right)
                        left = 1;
                        % RR = min(T, right + flankHalf);
                        % Rmax = max(sig(right:RR));
                        % cenMax = max(sig(left:right));
                        % prominence(ch,k,tr) = cenMax / max(Rmax, eps);

                    elseif ~isempty(left)
                        right = T; 
                        % LL = max(1, left - flankHalf);
                        % Lmax = max(sig(LL:left));
                        % cenMax = max(sig(left:right));
                        % prominence(ch,k,tr) = cenMax / max(Lmax, eps);

                    

                    end
                    if ~isempty(right) && ~isempty(left)
                        prominence(ch,k,tr) = local_prominence_snr( ...
                                    sig, left, right, flankHalf, 'pos', ...
                                    chGlobalMAD(ch), madFloorFrac);
                        % prominence(ch,k,tr) = cenMax / max(Rmax, Lmax, eps,...
                        %                         'omitnan');
                    else
                        prominence(ch,k,tr) = NaN; % no context available
                    end
                % --------- NEGATIVE spike ---------
                elseif sig(k) < lower
                    % first index to the RIGHT that rises back above 'lower'
                    right = find(sig(k:end) > lower, 1, 'first');
                    if ~isempty(right), right = k + right - 1; end

                    % last index to the LEFT that is above 'lower'
                    left  = find(sig(1:k) > lower, 1, 'last');

                    if ~isempty(right) && ~isempty(left)
                        % LL = max(1, left - flankHalf);
                        % Lmin = min(sig(LL:left));
                        % RR = min(T, right + flankHalf);
                        % Rmin = min(sig(right:RR));
                        % cenMin = min(sig(left:right));
                        % % use absolute values so ratio > 1 means "deeper" than flanks
                        % prominence(ch,k,tr) = abs(cenMin) / max([abs(Lmin), abs(Rmin), eps]);

                    elseif ~isempty(right)
                        left = 1;
                        % RR = min(T, right + flankHalf);
                        % Rmin = min(sig(right:RR));
                        % cenMin = min(sig(left:right));
                        % prominence(ch,k,tr) = abs(cenMin) / max(abs(Rmin), eps);

                    elseif ~isempty(left)
                        right = T; 
                        % LL = max(1, left - flankHalf);
                        % Lmin = min(sig(LL:left));
                        % cenMin = min(sig(left:right));
                        % prominence(ch,k,tr) = abs(cenMin) / max(abs(Lmin), eps);

                    
                    end
                    if ~isempty(right) && ~isempty(left)
                        prominence(ch,k,tr) = local_prominence_snr( ...
                                    sig, left, right, flankHalf, 'neg', ...
                                    chGlobalMAD(ch), madFloorFrac);
                        % prominence(ch,k,tr) = cenMax / max(Rmax, Lmax, eps,...
                        %                         'omitnan');
                    else
                        prominence(ch,k,tr) = NaN; % no context available
                    end
                end
            end
        end
    end
end






% If markWholeWindow==true, expand each detection to cover the full window samples around the detection center.
if ~markWholeWindow
    flags = double(detections);
else
    flags = zeros(C,T,N,'double');
    % Loop over channels and trials to expand detections (vectorized expansion is trickier but this is fine)
    for ch = 1:C
        for tr = 1:N
            idx = find(detections(ch,:,tr));
            for k = 1:numel(idx)
                t0 = idx(k);
                tstart = max(1, t0 - left);
                tstop  = min(T, t0 + right);
                flags(ch, tstart:tstop, tr) = 1;
            end
        end
    end
end

% If original input was 2D, squeeze outputs back to C x T
if origIs2D
    flags      = squeeze(flags);
    prominence = squeeze(prominence);
end

end
