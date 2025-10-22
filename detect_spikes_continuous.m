function [flags, prominence] = detect_spikes_continuous(data, threshold, ...
                    winSize, markWholeWindow, gammaDat)

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

[C,T,N] = size(data);

% Threshold broadcasting
if isscalar(threshold)
    thrMat = threshold; % scalar, MATLAB will broadcast
else
    threshold = threshold(:);
    if numel(threshold) ~= C
        error('If threshold is a vector it must have length C (channels).');
    end
    thrMat = repmat(reshape(threshold,[C 1 1]), 1, T, N); % C x T x N
end

% Moving max-min over time (dim=2)
half = floor((winSize-1)/2);
winLeft  = half;
winRight = winSize - 1 - winLeft; % same as half for odd win

rangeVals = movmax(data, [winLeft winRight], 2) - movmin(data, [winLeft winRight], 2);

% Detections
detections = rangeVals > thrMat;

% ---------- Prominence SNR ----------
prominence = zeros(size(detections), 'like', data);

% Per-channel global MAD (scaled)
chGlobalMAD = zeros(C,1);
for ch = 1:C
    x = reshape(gammaDat(ch,:,:), 1, []);
    chGlobalMAD(ch) = mad(x, 0); % scaled MAD
end
flankHalf     = 50;   % ~100-sample total flank window
madFloorFrac  = 0.3;  % floor as fraction of global MAD

for ch = 1:C
    % Channel-wise 10th/90th over all time×trials (works for 2D/3D)
    xall   = real(reshape(gammaDat(ch,:,:), 1, []));
    lims   = prctile(xall, [10 90]);
    lower  = lims(1); upper = lims(2);

    for tr = 1:N
        sig = squeeze(gammaDat(ch,:,tr));  % 1 x T

        for k = 1:T
            if ~detections(ch,k,tr), continue; end
            v = sig(k);
            if ~(v < lower || v > upper), continue; end

            if v > upper
                % positive: bracket by crossings below 'upper'
                rEdge = find(sig(k:end) < upper, 1, 'first');
                if ~isempty(rEdge), rEdge = k + rEdge - 1; end
                lEdge = find(sig(1:k)   < upper, 1, 'last');

                if isempty(lEdge) && ~isempty(rEdge), lEdge = 1; end
                if isempty(rEdge) && ~isempty(lEdge), rEdge = T; end

                if ~isempty(lEdge) && ~isempty(rEdge)
                    prominence(ch,k,tr) = local_prominence_snr( ...
                        sig, lEdge, rEdge, flankHalf, 'pos', ...
                        chGlobalMAD(ch), madFloorFrac);
                else
                    prominence(ch,k,tr) = NaN;
                end

            elseif v < lower
                % negative: bracket by crossings above 'lower'
                rEdge = find(sig(k:end) > lower, 1, 'first');
                if ~isempty(rEdge), rEdge = k + rEdge - 1; end
                lEdge = find(sig(1:k)   > lower, 1, 'last');

                if isempty(lEdge) && ~isempty(rEdge), lEdge = 1; end
                if isempty(rEdge) && ~isempty(lEdge), rEdge = T; end

                if ~isempty(lEdge) && ~isempty(rEdge)
                    prominence(ch,k,tr) = local_prominence_snr( ...
                        sig, lEdge, rEdge, flankHalf, 'neg', ...
                        chGlobalMAD(ch), madFloorFrac);
                else
                    prominence(ch,k,tr) = NaN;
                end
            end
        end
    end
end

% ---------- Expand detections to full window if requested ----------
if ~markWholeWindow
    flags = double(detections);
else
    flags = zeros(C,T,N,'double');
    for ch = 1:C
        for tr = 1:N
            idx = find(detections(ch,:,tr));
            for m = 1:numel(idx)
                t0 = idx(m);
                t1 = max(1, t0 - winLeft);
                t2 = min(T, t0 + winRight);
                flags(ch, t1:t2, tr) = 1;
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
