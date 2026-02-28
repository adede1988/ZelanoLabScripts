function chanDat = add_OB_gamma_features_to_behDat(chanDat)
% Adds OB gamma + burst timing/state features into chanDat.behDat.
% Assumes:
%   chanDat.fooof.spectra_flat_log10 : [nBreaths x 5 x nF]
%   chanDat.tf.powZ                  : [nBreaths x 50 x nF]
%   chanDat.tf.phase                 : [nBreaths x 50 x nF]
%   chanDat.tf.frex                  : [nF x 1]
%   chanDat.use                      : [nBreaths x 1] (1=QC pass)
%   chanDat.behDat                   : table with a condition field (see getCondVector below)

% -------------------- basics --------------------
bd   = chanDat.behDat;
useV = chanDat.use(:) == 1;
nB   = numel(useV);

frex = double(chanDat.tf.frex(:));
flat = double(chanDat.fooof.spectra_flat_log10); % [nB x 5 x nF]
powZ = double(chanDat.tf.powZ);                  % [nB x 50 x nF]
phTF = double(chanDat.tf.phase);                 % [nB x 50 x nF]
fooofGam = double(chanDat.fooof.gamma_peaks); 
nEpochs = 5;
nBins   = size(powZ,2);

% Ensure behDat has 'use' column (export requirement)
bd.use = single(useV);

% Condition vector for history reset
condVec = getCondVector(bd, nB, chanDat);

% Epoch naming -> column name prefixes
epPrefix = ["inRise","inFall","exRise","exFall","pause"];

% -------------------- (A) Gaussian fits per breath x epoch (25–60 Hz) --------------------
fitBand = [25 60];
sigmaBounds = [1 10];

maskFit = frex >= fitBand(1) & frex <= fitBand(2);
fFit = frex(maskFit);

freqFit     = nan(nB,nEpochs,'single');
fitR2       = nan(nB,nEpochs,'single');
promWide    = nan(nB,nEpochs,'single');
promLocal   = nan(nB,nEpochs,'single');
fooofBased  = zeros(nB,nEpochs,'single'); 
maskWide = frex >= 15 & frex <= 90;
maskLocalBase = frex >= 25 & frex <= 58;

for b = 1:nB
    for e = 1:nEpochs
        yAll = squeeze(flat(b,e,:));        % [nF x 1]
        yFit = yAll(maskFit);              % [nFit x 1]

        if all(~isfinite(yFit))
            continue
        end
        
        % Fit 1 Gaussian + offset on (fFit, yFit)
        out = fit_gauss1_offset_bounded(fFit, yFit, fitBand, sigmaBounds);
        if ~isnan(fooofGam(b,e))
            out.mu = fooofGam(b,e); 
            fooofBased(b,e) = 1; 
        end
        freqFit(b,e) = single(out.mu);
        fitR2(b,e)   = single(out.R2);

        if isfinite(out.mu)
            % value at fitted peak (use nearest frex bin in the full spectrum)
            [~, iMu] = min(abs(frex - out.mu));
            yMu = yAll(iMu);

            % promWide: baseline median 15–90
            baseW = median(yAll(maskWide), 'omitnan');
            promWide(b,e) = single(yMu - baseW);

            % promLocal: baseline median 25–58 excluding ±5 Hz around mu
            excl = abs(frex - out.mu) <= 5;
            baseMask = maskLocalBase & ~excl;
            baseL = median(yAll(baseMask), 'omitnan');
            promLocal(b,e) = single(yMu - baseL);
        end
    end
end

% Write epoch-wise columns
for e = 1:nEpochs
    bd.(epPrefix(e) + "Freq")     = freqFit(:,e);
    bd.(epPrefix(e) + "FitR2")    = fitR2(:,e);

    % Prominence outputs: both versions, plus keep the simple name as "wide"
    bd.(epPrefix(e) + "Prom")         = promWide(:,e);
    bd.(epPrefix(e) + "PromWide")     = promWide(:,e);
    bd.(epPrefix(e) + "PromLocal")    = promLocal(:,e);
    bd.(epPrefix(e) + "fooof")        = fooofBased(:,e);
    bd.allFooof                       = sum(fooofBased, 2, 'omitnan')>1; 
end

% Derived frequency summaries (per breath)
bd.minFreq    = single(min(freqFit, [], 2, 'omitnan'));
bd.maxFreq    = single(max(freqFit, [], 2, 'omitnan'));
bd.freqRange  = single(bd.maxFreq - bd.minFreq);

% Prom modulation index across epochs (from promWide)
promModIndex = nan(nB,1);
for b = 1:nB
    p = double(promWide(b,:));
    p = p(isfinite(p));
    if numel(p) >= 2
        mx = max(p); mn = min(p);
        den = (mx + mn);
        if den ~= 0
            promModIndex(b) = (mx - mn) / den;
        end
    end
end
bd.promModIndex = single(promModIndex);

% -------------------- (B) Participant center gamma freq + phi_burst features --------------------
% Participant center frequency from epoch fits, QC-only
centerHz = median(double(freqFit(useV,:)), 'all', 'omitnan');
if ~isfinite(centerHz)
    centerHz = mean(frex(frex>=25 & frex<=58), 'omitnan');
end

bandLo = prctile(freqFit(:), 25);
bandHi = prctile(freqFit(:), 75);
bandMask = frex >= bandLo & frex <= bandHi;
fBand = frex(bandMask);

% Max over freq within band at each time bin -> [nB x 50]
tmp = powZ(:,:,bandMask);
tmp(~isfinite(tmp)) = -Inf;
[maxGamma_t, argFi] = max(tmp, [], 3);         % max value and argmax (index within band)
maxGamma_t(maxGamma_t==-Inf) = NaN;

% Burst time index and value: max across time
[maxGamma_val, phi_idx] = max(maxGamma_t, [], 2); % [nB x 1] + [nB x 1]
phi_idx(~isfinite(maxGamma_val)) = NaN;

% Burst freq at that time (argmax freq within band at phi_idx)
phi_freq = nan(nB,1);
for b = 1:nB
    ti = phi_idx(b);
    if isfinite(ti)
        fi = argFi(b, ti);
        if isfinite(fi) && fi >= 1 && fi <= numel(fBand)
            phi_freq(b) = fBand(fi);
        end
    end
end

% Map idx (1..50) -> angle (-pi..pi) using bins-1 so that 1=-pi and 50=+pi
phi_ang = nan(nB,1);
for b = 1:nB
    ti = phi_idx(b);
    if isfinite(ti)
        phi_ang(b) = ((double(ti)-1) / (nBins-1)) * (2*pi) - pi;
    end
end

bd.gamma_centerHz = ones(nB, 1) * single(centerHz);
bd.gamma_bandLoHz = ones(nB, 1) * single(bandLo);
bd.gamma_bandHiHz = ones(nB, 1) * single(bandHi);

bd.phi_burst_idx  = single(phi_idx);
bd.phi_burst_freq = single(phi_freq);
bd.phi_burst      = single(phi_ang);
bd.cos_phi_burst  = single(cos(phi_ang));
bd.sin_phi_burst  = single(sin(phi_ang));

% prom_burst_f: peakVal - median across frequency at that timepoint
prom_burst_f = nan(nB,1);
prom_burst_t = nan(nB,1);

for b = 1:nB
    ti = phi_idx(b);
    if ~isfinite(ti) || ~isfinite(maxGamma_val(b)), continue; end

    baseF = median(squeeze(powZ(b,ti,:)), 'omitnan');
    prom_burst_f(b) = maxGamma_val(b) - baseF;

    baseT = median(maxGamma_t(b,:), 'omitnan');
    prom_burst_t(b) = maxGamma_val(b) - baseT;
end

bd.prom_burst_f = single(prom_burst_f);
bd.prom_burst_t = single(prom_burst_t);

% burstPeakWidthBins: threshold = 0.8 * median(prom_burst_t) across QC breaths
thr = 0.8 * median(prom_burst_t(useV), 'omitnan');
wBins = nan(nB,1);
if isfinite(thr)
    for b = 1:nB
        v = maxGamma_t(b,:);
        if all(~isfinite(v)), continue; end
        v0 = v - median(v, 'omitnan');        % breath-specific baseline (same as prom_burst_t definition)
        wBins(b) = sum(v0 >= thr);
    end
end
bd.burstPeakWidthBins = single(wBins);

% -------------------- (C) Rolling state features: prior 10 QC breaths within condition --------------------
phi_mean = nan(nB,1);
phi_err  = nan(nB,1);
R_win    = nan(nB,1);

mPromF   = nan(nB,1);
mPromT   = nan(nB,1);

for b = 1:nB
    c = condVec(b);

    % prior breaths in same condition and QC passed
    prior = find(useV & condVec==c & ( (1:nB)' < b ));
    if isempty(prior)
        continue
    end
    if numel(prior) > 10
        prior = prior(end-9:end);
    end

    ph = phi_ang(prior);
    ph = ph(isfinite(ph));
    if ~isempty(ph)
        z = mean(exp(1i*ph), 'omitnan');
        phi_mean(b) = angle(z);
        R_win(b)    = abs(z);
    end

    if isfinite(phi_ang(b)) && isfinite(phi_mean(b))
        phi_err(b) = angle(exp(1i*(phi_ang(b) - phi_mean(b))));
    end

    mPromF(b) = mean(prom_burst_f(prior), 'omitnan');
    mPromT(b) = mean(prom_burst_t(prior), 'omitnan');
end

bd.phi_mean          = single(phi_mean);
bd.phi_err           = single(phi_err);
bd.R_window          = single(R_win);
bd.mean_prom_burst_f = single(mPromF);
bd.mean_prom_burst_t = single(mPromT);

% -------------------- (D) Theta/alpha features at phi_burst_idx --------------------
% Theta: 2–10 Hz ; Alpha: 8.5–14 Hz
thetaBand = [2 10];
alphaBand = [8.5 14];
lowBaseBand = [1 25]; % for frequency-prom baselines

maskLowBase = frex>=lowBaseBand(1) & frex<=lowBaseBand(2);

theta_f = nan(nB,1); theta_f_prom = nan(nB,1); theta_t_prom = nan(nB,1); theta_phase = nan(nB,1);
alpha_f = nan(nB,1); alpha_f_prom = nan(nB,1); alpha_t_prom = nan(nB,1); alpha_phase = nan(nB,1);

maskTheta = frex>=thetaBand(1) & frex<=thetaBand(2);
maskAlpha = frex>=alphaBand(1) & frex<=alphaBand(2);

iTheta = find(maskTheta);
iAlpha = find(maskAlpha);

for b = 1:nB
    ti = phi_idx(b);
    if ~isfinite(ti), continue; end

    % ---- theta ----
    v = squeeze(powZ(b,ti,iTheta));
    v(~isfinite(v)) = -Inf;
    [mx, imx] = max(v);
    if mx ~= -Inf
        fi = iTheta(imx);
        theta_f(b) = frex(fi);
        theta_f_prom(b) = mx - median(squeeze(powZ(b,ti,maskLowBase)), 'omitnan');
        theta_t_prom(b) = mx - median(squeeze(powZ(b,:,fi)), 'omitnan');
        theta_phase(b)  = phTF(b,ti,fi);
    end

    % ---- alpha ----
    v = squeeze(powZ(b,ti,iAlpha));
    v(~isfinite(v)) = -Inf;
    [mx, imx] = max(v);
    if mx ~= -Inf
        fi = iAlpha(imx);
        alpha_f(b) = frex(fi);
        alpha_f_prom(b) = mx - median(squeeze(powZ(b,ti,maskLowBase)), 'omitnan');
        alpha_t_prom(b) = mx - median(squeeze(powZ(b,:,fi)), 'omitnan');
        alpha_phase(b)  = phTF(b,ti,fi);
    end
end

bd.burst_theta_f      = single(theta_f);
bd.burst_theta_f_prom = single(theta_f_prom);
bd.burst_theta_t_prom = single(theta_t_prom);
bd.burst_theta_phase  = single(theta_phase);

bd.burst_alpha_f      = single(alpha_f);
bd.burst_alpha_f_prom = single(alpha_f_prom);
bd.burst_alpha_t_prom = single(alpha_t_prom);
bd.burst_alpha_phase  = single(alpha_phase);



% -------------------- write back --------------------
chanDat.behDat = bd;

end

% ======================================================================
% Helper: condition vector (prefers behDat.condition; else best-effort)
% ======================================================================
function condVec = getCondVector(bd, nB, chanDat)
condVec = ones(nB,1);
try

    x = bd.condition;
    condVec = double(x(:));
 
catch
end


end

% ======================================================================
% Helper: bounded 1-Gaussian + offset fit in log10 space with fminsearch
% Model: yhat = A * exp(-(f-mu)^2/(2*sigma^2)) + C
% Constraints: mu in [fmin,fmax], sigma in [smin,smax]
% Allows A negative; C free
% Returns R2 on the fit points.
% ======================================================================
function out = fit_gauss1_offset_bounded(f, y, muBounds, sigBounds)

m = isfinite(f) & isfinite(y);
f = double(f(m));
y = double(y(m));

if numel(f) < 8
    out.mu = NaN; out.sigma = NaN; out.A = NaN; out.C = NaN; out.R2 = NaN;
    return
end

fmin = muBounds(1); fmax = muBounds(2);
smin = sigBounds(1); smax = sigBounds(2);

% init guesses
[~, im] = max(y);
mu0 = f(im);
C0  = median(y, 'omitnan');
A0  = y(im) - C0;
sig0 = 3;

% unconstrained parameters using sigmoids for bounds
p0 = [log(max(A0, eps)), ...
    inv_sigmoid01((mu0 - fmin) / (fmax - fmin)), ...
    inv_sigmoid01((sig0 - smin) / (smax - smin)), C0];

obj = @(p) sse_obj(p, f, y, fmin, fmax, smin, smax);

p = fminsearch(obj, p0, optimset('Display','off','MaxIter',600,'MaxFunEvals',2000));

[A, mu, sigma, C] = unpack_p(p, fmin, fmax, smin, smax);
yhat = A * exp(-0.5*((f - mu)./sigma).^2) + C;

SSE = sum((y - yhat).^2);
SST = sum((y - mean(y)).^2);

if SST > 0
    R2 = 1 - SSE/SST;
else
    R2 = NaN;
end

out.A = A; out.mu = mu; out.sigma = sigma; out.C = C; out.R2 = R2;

end

function J = sse_obj(p, f, y, fmin, fmax, smin, smax)
[A, mu, sigma, C] = unpack_p(p, fmin, fmax, smin, smax);
yhat = A * exp(-0.5*((f - mu)./sigma).^2) + C;
e = y - yhat;
J = sum(e.^2);
if ~isfinite(J), J = 1e18; end
end

function [A, mu, sigma, C] = unpack_p(p, fmin, fmax, smin, smax)
A = exp(p(1));
mu = fmin + (fmax - fmin) * sigmoid01(p(2));
sigma = smin + (smax - smin) * sigmoid01(p(3));
C = p(4);
end

function y = sigmoid01(x)
y = 1 ./ (1 + exp(-x));
end

function x = inv_sigmoid01(y)
y = min(max(y, 1e-6), 1-1e-6);
x = log(y ./ (1-y));
end
