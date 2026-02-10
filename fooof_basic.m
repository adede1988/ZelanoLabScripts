function out = fooof_basic(freq, psd, varargin)
% FOOOF_BASIC  FOOOF-like parameterization of a power spectrum (supports knee).
%
%   out = fooof_basic(freq, psd, 'Name', value, ...)
%
% Inputs
%   freq : [N x 1] frequency vector in Hz (must be > 0 in the fit range)
%   psd  : [N x 1] power spectrum (linear power). NaN/Inf are ignored.
%
% Name-Value options (defaults in brackets)
%   'f_range'           : [fmin fmax] Hz subset to fit ([min(freq) max(freq)])
%   'aperiodic_mode'    : 'knee' or 'fixed'  ['knee']
%   'max_peaks'         : maximum number of peaks to fit [6]
%   'peak_thresh'       : minimum peak amplitude (SD of flattened spectrum) [2.0]
%   'peak_width_limits' : [fwhm_min fwhm_max] in Hz for Gaussian width [1 12]
%   'smooth_w'          : odd integer window (points) to smooth flattened spec [5]
%   'max_clip_iters'    : sigma-clipping iterations for aperiodic fit [3]
%   'verbose'           : true/false [false]
%
% Outputs (struct)
%   out.ap.offset     : intercept b
%   out.ap.exponent   : exponent chi
%   out.ap.knee       : knee parameter k (only for 'knee' mode; else 0)
%   out.ap.knee_freq  : knee frequency in Hz (= k^(1/chi), NaN if undefined)
%   out.peaks         : [nPeaks x 3] = [center_Hz, amplitude_log10, fwhm_Hz]
%   out.fit_log10     : fields ap, peaks, full on the chosen fit range
%   out.flattened     : flattened spectrum (log10 power minus aperiodic) on range
%   out.meta          : options & bookkeeping
%
% Notes
% - Gaussians are fit/additive in log10 power (FOOOF style).
% - Aperiodic model can be refit once on the peak-removed spectrum.

% ---------- parse & validate ----------
p = inputParser;
p.addParameter('f_range', [min(freq(:)) max(freq(:))]);
p.addParameter('aperiodic_mode', 'knee');     % 'knee' or 'fixed'
p.addParameter('max_peaks', 6);
p.addParameter('peak_thresh', 2.0);
p.addParameter('peak_width_limits', [1 12]);  % Hz, FWHM
p.addParameter('smooth_w', 5);
p.addParameter('max_clip_iters', 3);
p.addParameter('verbose', false);
p.parse(varargin{:});
opt = p.Results;

freq = freq(:);
psd  = psd(:);
if numel(freq) ~= numel(psd), error('freq and psd must have same length'); end
if any(~isfinite(freq)), error('freq must be finite'); end

% restrict to fitting range and valid samples
in_rng = freq >= opt.f_range(1) & freq <= opt.f_range(2) & freq > 0 & isfinite(psd);
f  = freq(in_rng);
P  = psd(in_rng);
if numel(f) < 16, error('Not enough valid points in the specified f_range.'); end

% ---------- log transforms ----------
x = log10(f);              % predictor for fixed model
y = log10(P);              % response
% figure
% plot(x, y)
% ---------- initial aperiodic fit (robust, sigma-clipped) ----------
[ap_params, ap_fit] = fit_aperiodic(f, x, y, opt.aperiodic_mode, opt.max_clip_iters);

% flatten for peak finding
flat = y - ap_fit;

% keep a smoothed copy for *detection* (fit uses unsmoothed)
w = opt.smooth_w; if mod(w,2)==0, w = w+1; end
flat_s = movmean(flat, max(1, w));

% robust noise estimate (SD of flattened spectrum)
noise_sd = 1.4826 * mad(flat, 1);

% ---------- peak detection & fitting loop ----------
fwhm_min = opt.peak_width_limits(1);
fwhm_max = opt.peak_width_limits(2);
sig_min  = fwhm_min / (2*sqrt(2*log(2)));
sig_max  = fwhm_max / (2*sqrt(2*log(2)));

peaks = [];                 % rows = [mu_Hz, amp_log10, fwhm_Hz]
resid = flat;               % residual (flattened) to deflate as peaks are fit

nfit = 0;
while nfit < opt.max_peaks
    [iPeak, val] = local_max_idx(flat_s);
    if isempty(iPeak), break; end
    [valMax, ord] = max(val);
    i0 = iPeak(ord);
    if valMax < opt.peak_thresh * noise_sd
        break;
    end

    mu0   = f(i0);
    A0    = max(resid(i0), 0);
    sig0  = (sig_min + sig_max) / 2;

    % window around candidate (±3*sig_max)
    win = (f >= (mu0 - 3*sig_max)) & (f <= (mu0 + 3*sig_max));
    if nnz(win) < 7, flat_s(i0) = -Inf; continue; end

    fw   = f(win);
    rw   = resid(win);

    p0 = [max(A0, eps), mu0, max(sig0, eps)];
    cost = @(p) peak_cost(p, fw, rw, mu0, sig_min, sig_max);

    p_opt = fminsearch(cost, p0, optimset('Display','off'));
    A = max(p_opt(1), 0);
    mu = p_opt(2);
    sg = max(min(p_opt(3), sig_max), sig_min);
    fwhm = 2*sqrt(2*log(2)) * sg;

    if ~all(isfinite([A mu sg])) || A < opt.peak_thresh * 0.5 * noise_sd
        flat_s(i0) = -Inf;
        continue
    end

    g = A * exp(-0.5 * ((f - mu)./sg).^2);
    resid  = resid - g;
    flat_s = movmean(resid, max(1, w));

    peaks = [peaks; mu, A, fwhm]; %#ok<AGROW>
    nfit  = nfit + 1;

    if opt.verbose
        fprintf('Peak %d: mu=%.2f Hz, amp=%.3f log10, fwhm=%.2f Hz\n', nfit, mu, A, fwhm);
    end
end

% ----- build Gaussian design & initial peak sum on the range -----
K = size(peaks,1);
peak_fit = zeros(size(f));
if K > 0
    G = zeros(numel(f), K);
    for k = 1:K
        mu   = peaks(k,1);
        fwhm = peaks(k,3);
        sg   = fwhm / (2*sqrt(2*log(2)));
        G(:,k) = exp(-0.5*((f - mu)./sg).^2);
    end
    A = peaks(:,2);
    peak_fit = G * A;
end

% ----- refit aperiodic on peak-removed spectrum -----
y_nopeaks = y - peak_fit;
[ap_params2, ap_fit2] = fit_aperiodic(f, x, y_nopeaks, opt.aperiodic_mode, opt.max_clip_iters);

% ----- re-estimate peak amplitudes against the new baseline -----
r2 = y - ap_fit2;                      % flattened wrt updated aperiodic
if K > 0
    if exist('lsqnonneg','file')
        A2 = lsqnonneg(G, r2);
    else
        A2 = max(0, G \ r2);
    end
    peak_fit2 = G * A2;
    peaks(:,2) = A2;                   % update amplitudes
else
    peak_fit2 = zeros(size(f));
end
full_fit2 = ap_fit2 + peak_fit2;

% ---------- package outputs ----------
% aperiodic params (mode-dependent)
switch lower(opt.aperiodic_mode)
    case 'fixed'
        b   = ap_params2(1);
        chi = ap_params2(2);
        k   = 0;
        knee_freq = NaN;
    case 'knee'
        b   = ap_params2(1);
        k   = 10.^ap_params2(2);       % params store log10(k)
        chi = ap_params2(3);
        if chi > 0 && isfinite(k)
            knee_freq = k^(1/chi);
        else
            knee_freq = NaN;
        end
    otherwise
        error('Unknown aperiodic_mode: %s', opt.aperiodic_mode);
end

out.ap.offset     = b;
out.ap.exponent   = chi;
out.ap.knee       = k;
out.ap.knee_freq  = knee_freq;

out.peaks             = peaks;  % [center_Hz, amplitude_log10, fwhm_Hz]
out.fit_log10.ap      = ap_fit2;
out.fit_log10.peaks   = peak_fit2;
out.fit_log10.full    = full_fit2;
out.flattened         = r2;

out.meta.f_range           = opt.f_range;
out.meta.aperiodic_mode    = opt.aperiodic_mode;
out.meta.max_peaks         = opt.max_peaks;
out.meta.peak_thresh       = opt.peak_thresh;
out.meta.peak_width_limits = opt.peak_width_limits;
out.meta.noise_sd          = noise_sd;

% ----- also return full-length reconstructions -----
out.full.freq = freq(:);
switch lower(opt.aperiodic_mode)
    case 'fixed'
        ap_full = b - chi * log10(max(out.full.freq, eps));
    case 'knee'
        ap_full = b - log10(10.^ap_params2(2) + max(out.full.freq, eps).^chi);
end
out.full.ap_log10     = ap_full;
out.full.peaks_log10  = zeros(size(out.full.freq));
for kix = 1:size(peaks,1)
    A = peaks(kix,2); mu = peaks(kix,1); fwhm = peaks(kix,3);
    sg = fwhm / (2*sqrt(2*log(2)));
    out.full.peaks_log10 = out.full.peaks_log10 + A * exp(-0.5*((out.full.freq - mu)./sg).^2);
end
out.full.model_log10 = out.full.ap_log10 + out.full.peaks_log10;

end

% ======== helper: robust aperiodic fit (fixed or knee) ========
function [params, ap_fit] = fit_aperiodic(f, x, y, mode, max_clip_iters)
% Returns:
%   params:
%     fixed: [b, chi]
%     knee : [b, log10k, chi]
%   ap_fit: model prediction at the provided f/x

mode = lower(mode);
mask = isfinite(y);
% initial guess from linear model
co = polyfit(x(mask), y(mask), 1);    % y ≈ m*x + b
m0 = co(1); b0 = co(2);
chi0 = -m0;
log10k0 = log10(max( (min(f)+eps)^max(chi0,0), 1e-6 )); % conservative tiny knee

switch mode
    case 'fixed'
        p = [b0, max(chi0, 0)];       % [b, chi]
    case 'knee'
        p = [b0, log10k0, max(chi0, 0.1)];   % [b, log10k, chi]
    otherwise
        error('Unknown aperiodic_mode: %s', mode);
end

for it = 1:max_clip_iters
    switch mode
        case 'fixed'
            yhat = p(1) - p(2)*x;
        case 'knee'
            k = 10.^p(2);
            yhat = p(1) - log10(k + f.^p(3));
    end
    r = y - yhat;
    s = nanstd(r(mask));
    if ~isfinite(s) || s==0, break; end
    newmask = mask & (abs(r) <= 2.5*s);
    if isequal(newmask, mask), break; end
    mask = newmask;

    % refit on the clipped mask with simple optimizer
    obj = @(pp) ap_obj(pp, f(mask), x(mask), y(mask), mode);
    p = fminsearch(obj, p, optimset('Display','off'));
end

% final refine on current mask
obj = @(pp) ap_obj(pp, f(mask), x(mask), y(mask), mode);
p   = fminsearch(obj, p, optimset('Display','off'));

switch mode
    case 'fixed'
        params = [p(1), p(2)];
        ap_fit = p(1) - p(2)*x;
    case 'knee'
        params = [p(1), p(2), p(3)];  % [b, log10k, chi]
        k = 10.^p(2);
        ap_fit = p(1) - log10(k + f.^p(3));
end
end

function J = ap_obj(p, f, x, y, mode)
% objective with soft bounds/penalties for stability
mode = lower(mode);
switch mode
    case 'fixed'
        b = p(1); chi = p(2);
        yhat = b - chi*x;
        pen = 0;
        % gentle bounds for chi
        if ~isfinite(chi) || chi < 0 || chi > 6
            pen = pen + 1e6*(~isfinite(chi) + max(0, -chi) + max(0, chi-6)).^2;
        end
    case 'knee'
        b = p(1); log10k = p(2); chi = p(3);
        k = 10.^log10k;
        yhat = b - log10(k + f.^chi);
        pen = 0;
        if ~isfinite(chi) || chi < 0.1 || chi > 6
            pen = pen + 1e6*(~isfinite(chi) + max(0, 0.1-chi) + max(0, chi-6)).^2;
        end
        if ~isfinite(log10k) || log10k < -12 || log10k > 12  % 1e-12 .. 1e12 in k-units
            pen = pen + 1e6*(~isfinite(log10k) + max(0, -12-log10k) + max(0, log10k-12)).^2;
        end
end
e  = y - yhat;
J  = sum(e.^2) + pen;
end

% ===== helper: local maxima indices (simple, no toolbox) =====
function [idx, vals] = local_max_idx(v)
n = numel(v);
if n < 3, idx = []; vals = []; return; end
is_peak = false(n,1);
for i = 2:n-1
    if v(i) > v(i-1) && v(i) >= v(i+1)
        is_peak(i) = true;
    end
end
idx  = find(is_peak);
vals = v(idx);
end

% ===== helper: Gaussian fit cost with soft bounds =====
function J = peak_cost(p, f, r, mu_hint, sig_min, sig_max)
A  = p(1); mu = p(2); sg = p(3);
g  = A * exp(-0.5 * ((f - mu)./sg).^2);
e  = r - g;
J  = sum(e.^2);

if ~isfinite(A) || ~isfinite(mu) || ~isfinite(sg) || sg <= 0
    J = J + 1e9; return;
end
if A < 0
    J = J + 1e6 * (abs(A));
end
if sg < sig_min
    J = J + 1e6 * (sig_min - sg)^2;
elseif sg > sig_max
    J = J + 1e6 * (sg - sig_max)^2;
end
if abs(mu - mu_hint) > 2*sig_max
    J = J + 1e4 * (abs(mu - mu_hint) - 2*sig_max)^2;
end
end
