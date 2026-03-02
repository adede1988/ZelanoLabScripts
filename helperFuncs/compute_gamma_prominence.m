
function [promEpoch, promPooled] = compute_gamma_prominence(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
% Returns median prominence per epoch and pooled median across all epochs
[promByEpoch, promAll] = compute_gamma_prominence_full(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak);
nEpochs = numel(promByEpoch);

promEpoch = nan(1,nEpochs);
for e=1:nEpochs
    promEpoch(e) = median(promByEpoch{e}, 'omitnan');
end
promPooled = median(promAll, 'omitnan');
end