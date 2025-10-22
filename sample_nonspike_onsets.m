function nonSpikeOnsets = sample_nonspike_onsets(spikeOnsets, T, pre, post, seed)
% SAMPLE_NONSPIKE_ONSETS
%   spikeOnsets : M x 1 (sample indices of spike onsets, 1-based)
%   T           : total number of time samples in the recording
%   pre, post   : window half-widths you plan to extract around each onset
%   seed        : (optional) RNG seed for reproducibility
%
%   Returns:
%     nonSpikeOnsets : M x 1 indices, each a valid center for a ±[pre,post] window
%
%   Notes:
%     - Ensures centers are ≥ pre from the start and ≤ T-post from the end
%     - Excludes any index whose ±[pre,post] window would overlap a spike window
%     - Samples without replacement if enough candidates exist; otherwise with replacement

if nargin >= 5 && ~isempty(seed), rng(seed); end

spikeOnsets = spikeOnsets(:);
M = numel(spikeOnsets);

% 1) Start from all valid centers that can support a full window
good = true(1, T);
good(1:pre) = false;
good(T-post+1:T) = false;

% 2) Exclude windows around each spike onset (so baselines are spike-free)
for s = spikeOnsets.'
    i1 = max(1, s - pre);
    i2 = min(T, s + post);
    good(i1:i2) = false;
end

candidates = find(good);
numCand = numel(candidates);

if numCand == 0
    error('No valid non-spike centers remain with the given pre/post exclusions.');
end

% 3) Sample M centers
if numCand >= M
    % without replacement
    pick = randperm(numCand, M);
    nonSpikeOnsets = sort(candidates(pick)).';
else
    % with replacement (duplicates possible if the pool is small)
    idx = randi(numCand, M, 1);
    nonSpikeOnsets = sort(candidates(idx));
end
end
