function ridgeInfo = v2_ridges(T, C)
% V2_RIDGES  Primary + secondary ridges on the z-scored FASLT TFR (spec 4).
%   ridgeInfo = v2_ridges(T, C)   % T from v2_tfr
%
%   Primary: DP ridge (tfridge; penalty C.ridge.penalty) within C.ridge.band on the z-scored
%   TFR. Peel: at each time along the ridge, fit a Gaussian to the primary peak using its FWHM
%   extent along frequency and subtract it from the TFR; re-extract the secondary ridge from the
%   residual. Noisy/invalid trials stay NaN so matrices stay aligned to behDat rows.
%
%   (tfridge is the MATLAB-native DP ridge tracker sanctioned by the spec; the recommended
%   Frequency_ridge_tracking repo is Python-only, MODA ecurve is the MATLAB fallback.)
%
%   ridgeInfo.primaryRidge.f/.p  [N x nCore]   (Hz / z-power along ridge)
%   ridgeInfo.secondaryRidge.f/.p[N x nCore]
%   ridgeInfo.tMs freqs band penalty good nTrials

    freqs = T.freqs(:)'; nCore = T.nCore; N = T.N;
    bandIdx = find(freqs >= C.ridge.band(1) & freqs <= C.ridge.band(2));
    fBand = freqs(bandIdx);
    pen = C.ridge.penalty; minPk = C.ridge.minPeakZ;

    pf = nan(N,nCore); pp = nan(N,nCore); sf = nan(N,nCore); sp = nan(N,nCore);

    for i = find(T.good)'
        Z  = T.zTFR(:,:,i);              % [Nf x nCore] full freq (signed z)
        if all(~isfinite(Z(:))), continue; end
        Zb = Z(bandIdx,:);
        % --- primary ridge (search on nonneg energy; record true z along ridge) ---
        [fr, ir] = tfridge(max(Zb,0), fBand, pen, 'NumRidges', 1);
        fr = fr(:)'; ir = ir(:)';
        pr = Zb(sub2ind(size(Zb), ir, 1:nCore));
        pf(i,:) = fr; pp(i,:) = pr;

        % --- Gaussian FWHM peel of the primary, per time point ---
        R = Z;
        if C.ridge.fwhmPeel
            for t = 1:nCore
                fullc = bandIdx(ir(t));           % full-freq index of the ridge at time t
                pk = Z(fullc, t);
                if ~(pk > minPk), continue; end   % skip noise-floor times
                s = fwhmSigma(Z(:,t), fullc, pk, freqs);
                G = pk .* exp(-((freqs - freqs(fullc)).^2) ./ (2*s^2));
                R(:,t) = Z(:,t) - G(:);
            end
        end
        % --- secondary ridge from the residual ---
        Rb = R(bandIdx,:);
        [fr2, ir2] = tfridge(max(Rb,0), fBand, pen, 'NumRidges', 1);
        fr2 = fr2(:)'; ir2 = ir2(:)';
        pr2 = Rb(sub2ind(size(Rb), ir2, 1:nCore));
        sf(i,:) = fr2; sp(i,:) = pr2;
    end

    ridgeInfo = struct();
    ridgeInfo.primaryRidge   = struct('f', pf, 'p', pp);
    ridgeInfo.secondaryRidge = struct('f', sf, 'p', sp);
    ridgeInfo.tMs   = T.tMs;
    ridgeInfo.freqs = freqs;
    ridgeInfo.band  = C.ridge.band;
    ridgeInfo.penalty = pen;
    ridgeInfo.good  = T.good;
    ridgeInfo.nTrials = N;
    ridgeInfo.note  = 'primary/secondary ridge on z-scored FASLT TFR; tfridge + Gaussian-FWHM peel';
end

% ---- FWHM -> gaussian sigma (Hz), robust to edges / no-crossing ----
function s = fwhmSigma(col, ic, pk, freqs)
    half = pk/2;
    iu = ic; while iu < numel(col) && isfinite(col(iu+1)) && col(iu+1) >= half, iu = iu+1; end
    id = ic; while id > 1         && isfinite(col(id-1)) && col(id-1) >= half, id = id-1; end
    fwhmHz = freqs(iu) - freqs(id);
    if ~(fwhmHz > 0) || ~isfinite(fwhmHz), fwhmHz = 3; end   % default ~3 Hz bump
    s = fwhmHz / (2*sqrt(2*log(2)));
    if s < 0.5, s = 0.5; end
end
