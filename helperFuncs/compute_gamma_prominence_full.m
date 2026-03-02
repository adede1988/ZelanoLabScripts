

function [promByEpoch, promAll, promRaw] = compute_gamma_prominence_full(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak)
% Returns per-epoch prominence samples (cell), pooled sample vector, and raw per-breath×epoch matrix

Gdet = double(chanDat.fooof.gamma_peaks);        % NaN if no peak
Gmax = double(chanDat.fooof.gamma_peak_freq);    % fallback
flat = double(chanDat.fooof.spectra_flat_log10); % breaths x epochs x frex
frex = double(chanDat.tf.frex(:));

nBreaths = size(Gdet,1);
nEpochs  = size(Gdet,2);

% Choose peak frequency per breath/epoch: detected if possible else max
Guse = Gdet;
missing = ~isfinite(Guse);
Guse(missing) = Gmax(missing);

% Precompute masks for baseline band
maskBaseBand = frex>=baselineBandHz(1) & frex<=baselineBandHz(2);

promRaw = nan(nBreaths, nEpochs);
for b=1:nBreaths
    if ~useVec(b), continue; end
    for e=1:nEpochs
        fpk = Guse(b,e);
        if ~isfinite(fpk), continue; end
        [~, fiPk] = min(abs(frex - fpk));

        % baseline mask excludes +/- excludeHzAroundPeak around peak
        maskEx = maskBaseBand & ~(frex >= (fpk-excludeHzAroundPeak) & frex <= (fpk+excludeHzAroundPeak));

        base = median(flat(b,e,maskEx), 3, 'omitnan');
        if ~isfinite(base), continue; end

        promRaw(b,e) = flat(b,e,fiPk) - base;
    end
end

promByEpoch = cell(1,nEpochs);
promAll = [];
for e=1:nEpochs
    x = promRaw(useVec, e);
    x = x(isfinite(x));
    promByEpoch{e} = x(:);
    promAll = [promAll; x(:)];
end

end

