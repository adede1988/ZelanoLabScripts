function [Xlock, onsetSamples, onsetTrials] = lock_power_to_blinks(blinkMask, P, pre, post)
% LOCK_POWER_TO_BLINKS
%   Inputs:
%     blinkMask : T x N   (1 where blink present at that time; 0 otherwise)
%     P         : T x N x F  (power values)
%     pre       : samples before onset to include (default 500)
%     post      : samples after onset to include (default 500)
%
%   Outputs:
%     Xlock        : (pre+post+1) x B x F   (time x blinks x frequency)
%     onsetSamples : B x 1   onset sample indices (within trial)
%     onsetTrials  : B x 1   trial index for each blink

if nargin < 3 || isempty(pre),  pre  = 500; end
if nargin < 4 || isempty(post), post = 500; end

% --- basic checks ---
[T1, N1] = size(blinkMask);
[T2, N2, F] = size(P);
if T1 ~= T2 || N1 ~= N2
    error('blinkMask (T x N) and P (T x N x F) must share T and N.');
end
% logical-ize mask
blinkMask = blinkMask ~= 0;

% --- find valid onsets (0->1), keeping only those with full ±window ---
onsetSamples = [];
onsetTrials  = [];
for tr = 1:N1
    v = blinkMask(:,tr);
    on = find(diff([false; v]) == 1);       % indices where 0->1
    on = on(on > pre & on <= (T1 - post));  % keep only those with full window
    if ~isempty(on)
        onsetSamples = [onsetSamples; on(:)];
        onsetTrials  = [onsetTrials;  repmat(tr, numel(on), 1)];
    end
end

% --- build locked tensor: time x blinks x frequency ---
B = numel(onsetSamples);
W = pre + post + 1;
Xlock = zeros(W, B, F, 'like', P);

for b = 1:B
    t0 = onsetSamples(b);
    tr = onsetTrials(b);
    Xlock(:, b, :) = P(t0-pre : t0+post, tr, :);
end
end
