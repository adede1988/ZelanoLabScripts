function out = fooof_basic(freq, psd, varargin)
% FOOOF_BASIC  FOOOF-like parameterization of a power spectrum (supports knee).
% Iterative refitting version: alternates aperiodic + peak model updates.
%
%   out = fooof_basic(freq, psd, 'Name', value, ...)
%
% Inputs
%   freq : [N x 1] frequency vector in Hz
%   psd  : [N x 1] power spectrum (linear power). NaN/Inf are ignored.
%
% Name-Value options (defaults)
%   'f_range'           : [fmin fmax] Hz subset to fit ([min(freq) max(freq)])
%   'aperiodic_mode'    : 'knee' or 'fixed'  ['knee']
%   'max_peaks'         : maximum number of peaks to fit [6]
%   'peak_thresh'       : minimum peak amplitude (SD of flattened spectrum) [2.0]
%   'peak_width_limits' : [fwhm_min fwhm_max] in Hz for Gaussian width [1 12]
%   'smooth_w'          : odd integer window (points) to smooth flattened spec [5]
%   'max_clip_iters'    : sigma-clipping iterations for aperiodic fit [3]
%   'clip_sigma'        : sigma threshold for clipping [2.5]
%   'n_iter'            : OUTER iterations of (ap -> peaks -> ap -> peaks) [2]
%   'tol'               : convergence tolerance on aperiodic params [1e-3]
%   'peak_mask_mult'    : mask width for refitting ap (multiplier on sigma) [1.5]
%   'verbose'           : true/false [false]
%
% Outputs (struct)
%   out.ap.offset     : intercept b
%   out.ap.exponent   : exponent chi
%   out.ap.knee       : knee parameter k (only for 'knee'; else 0)
%   out.ap.knee_freq  : k^(1/chi), NaN if undefined
%   out.peaks         : [nPeaks x 3] = [center_Hz, amplitude_log10, fwhm_Hz]
%   out.fit_log10     : fields ap, peaks, full on the chosen fit range
%   out.flattened     : flattened spectrum (log10 power minus aperiodic) on range
%   out.meta          : options & bookkeeping
%   out.full          : full-length reconstructions (log10) on original freq grid
%
% Notes
% - Peaks are additive Gaussians in log10 power (FOOOF style).
% - Aperiodic model is refit on peak-removed spectrum; this is repeated n_iter times.

% ---------- parse & validate ----------
p = inputParser;
p.addParameter('f_range', [min(freq(:)) max(freq(:))]);
p.addParameter('aperiodic_mode', 'knee');     % 'knee' or 'fixed'
p.addParameter('max_peaks', 6);
p.addParameter('peak_thresh', 2.0);
p.addParameter('peak_width_limits', [1 12]);  % Hz, FWHM
p.addParameter('smooth_w', 5);
p.addParameter('max_clip_iters', 3);
p.addParameter('clip_sigma', 2.5);
p.addParameter('n_iter', 2);
p.addParameter('tol', 1e-3);
p.addParameter('peak_mask_mult', 1);
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

if numel(f) < 16
    error('Not enough valid points in the specified f_range.');
end

% guard power
P = double(P);
P(P <= 0) = eps;

% log transforms
x = log10(f);
y = log10(P);

% smoothing window (index space; for log-spaced frex this is ~log-frequency smoothing)
w = opt.smooth_w; if mod(w,2)==0, w = w+1; end
w = max(1, w);

% peak width bounds (Hz)
fwhm_min = opt.peak_width_limits(1);
fwhm_max = opt.peak_width_limits(2);
sig_min  = fwhm_min / (2*sqrt(2*log(2)));
sig_max  = fwhm_max / (2*sqrt(2*log(2)));

% ---------- iterative outer loop ----------
peaks = [];                 % [mu, amp, fwhm] on fit-range f
ap_params = [];             % mode-specific internal params
ap_fit = nan(size(f));      % aperiodic fit in log10 power on f
peak_fit = zeros(size(f));  % summed gaussians in log10 power on f

ap_prev = [];
converged = false;

for itOuter = 1:max(1, opt.n_iter)

    % 1) Build peak mask for aperiodic fitting (exclude regions around current peaks)
    peakMask = false(size(f));
    if ~isempty(peaks)
        for k = 1:size(peaks,1)
            mu   = peaks(k,1);
            fwhm = peaks(k,3);
            sg   = fwhm / (2*sqrt(2*log(2)));
            rad  = opt.peak_mask_mult * 3 * sg;   % ~±(3*sigma) scaled
            peakMask = peakMask | (f >= (mu-rad) & f <= (mu+rad));
        end
    end

    % 2) Fit aperiodic on (y - current peak_fit) using robust sigma clipping, excluding peakMask
    y_for_ap = y - peak_fit;
    [ap_params, ap_fit] = fit_aperiodic_masked(f, x, y_for_ap, opt.aperiodic_mode, ...
                                               opt.max_clip_iters, opt.clip_sigma, peakMask);

    % 3) Flatten and estimate noise
    flat = y - ap_fit;                           % flattened wrt updated aperiodic (not removing peaks yet)
    noise_sd = 1.4826 * mad(flat(isfinite(flat)), 1);
    if ~isfinite(noise_sd) || noise_sd <= 0
        noise_sd = std(flat(isfinite(flat)));
        if ~isfinite(noise_sd) || noise_sd <= 0, noise_sd = 1; end
    end

    % 4) Detect+fit peaks ITERATIVELY on residual flattened spectrum
    %    Canonical-ish: fit peaks on flattened, subtract, repeat.
    resid = flat;                                % what we deflate with Gaussians
    resid_s = movmean(resid, w);

    peaks_new = [];
    nfit = 0;

    while nfit < opt.max_peaks
        [idx, vals] = local_max_idx(resid_s);
        if isempty(idx), break; end

        % pick strongest peak
        [valMax, ord] = max(vals);
        i0 = idx(ord);

        % threshold in SD units
        if valMax < opt.peak_thresh * noise_sd
            break;
        end

        mu0  = f(i0);
        A0   = max(resid(i0), 0);               % amplitude in log10 units
        sg0  = (sig_min + sig_max)/2;

        % window around candidate
        win = (f >= (mu0 - 3*sig_max)) & (f <= (mu0 + 3*sig_max));
        if nnz(win) < 7
            resid_s(i0) = -Inf;
            continue;
        end

        fw = f(win);
        rw = resid(win);

        p0   = [max(A0, eps), mu0, max(sg0, eps)];
        cost = @(pp) peak_cost(pp, fw, rw, mu0, sig_min, sig_max);

        p_opt = fminsearch(cost, p0, optimset('Display','off'));

        A  = max(p_opt(1), 0);
        mu = p_opt(2);
        sg = max(min(p_opt(3), sig_max), sig_min);
        fwhm = 2*sqrt(2*log(2)) * sg;

        % reject degenerate
        if ~all(isfinite([A mu sg])) || A < (0.5*opt.peak_thresh*noise_sd)
            resid_s(i0) = -Inf;
            continue;
        end

        % subtract this peak from residual
        g = A * exp(-0.5 * ((f - mu)./sg).^2);
        resid   = resid - g;
        resid_s = movmean(resid, w);

        peaks_new = [peaks_new; mu, A, fwhm]; %#ok<AGROW>
        nfit = nfit + 1;

        if opt.verbose
            fprintf('[%d] Peak %d: mu=%.2f Hz, amp=%.3f log10, fwhm=%.2f Hz\n', ...
                itOuter, nfit, mu, A, fwhm);
        end
    end

    % 5) Build Gaussian design for these peaks and refit peak amplitudes against current aperiodic
    peaks = peaks_new;
    K = size(peaks,1);

    if K > 0
        G = zeros(numel(f), K);
        for k = 1:K
            mu   = peaks(k,1);
            fwhm = peaks(k,3);
            sg   = fwhm / (2*sqrt(2*log(2)));
            G(:,k) = exp(-0.5*((f - mu)./sg).^2);
        end

        % target residual wrt aperiodic
        r = y - ap_fit;

        % NNLS for nonnegative peak amps
        if exist('lsqnonneg','file') == 2
            A = lsqnonneg(G, r);
        else
            A = max(0, G \ r);
        end

        peak_fit = G * A;
        peaks(:,2) = A; % update amps
    else
        peak_fit = zeros(size(f));
    end

    % 6) Refit aperiodic on peak-removed spectrum (canonical step)
    y_nopeaks = y - peak_fit;

    % mask again using current peaks (more stable)
    peakMask = false(size(f));
    if ~isempty(peaks)
        for k = 1:size(peaks,1)
            mu   = peaks(k,1);
            fwhm = peaks(k,3);
            sg   = fwhm / (2*sqrt(2*log(2)));
            rad  = opt.peak_mask_mult * 3 * sg;
            peakMask = peakMask | (f >= (mu-rad) & f <= (mu+rad));
        end
    end

    [ap_params2, ap_fit2] = fit_aperiodic_masked(f, x, y_nopeaks, opt.aperiodic_mode, ...
                                                 opt.max_clip_iters, opt.clip_sigma, peakMask);

    % update aperiodic
    ap_fit = ap_fit2;
    ap_params = ap_params2;

    % 7) Re-estimate peak amplitudes one more time against updated aperiodic (canonical-ish)
    if K > 0
        r2 = y - ap_fit;
        if exist('lsqnonneg','file') == 2
            A2 = lsqnonneg(G, r2);
        else
            A2 = max(0, G \ r2);
        end
        peak_fit = G * A2;
        peaks(:,2) = A2;
    end

    % 8) Convergence check (aperiodic params stable)
    ap_now = ap_params(:);
    if ~isempty(ap_prev)
        if all(isfinite(ap_now)) && all(isfinite(ap_prev))
            if norm(ap_now - ap_prev) / max(1, norm(ap_prev)) < opt.tol
                converged = true;
                if opt.verbose
                    fprintf('Converged at outer iter %d\n', itOuter);
                end
                break;
            end
        end
    end
    ap_prev = ap_now;

end % outer loop

full_fit = ap_fit + peak_fit;
flat_out = y - ap_fit;               % flattened wrt final aperiodic (FOOOF-style)

% ---------- package outputs ----------
% decode aperiodic params to canonical outputs
switch lower(opt.aperiodic_mode)
    case 'fixed'
        b   = ap_params(1);
        chi = ap_params(2);
        k   = 0;
        knee_freq = NaN;
    case 'knee'
        b   = ap_params(1);
        k   = 10.^ap_params(2);      % stored as log10(k)
        chi = ap_params(3);
        if chi > 0 && isfinite(k)
            knee_freq = k^(1/chi);
        else
            knee_freq = NaN;
        end
    otherwise
        error('Unknown aperiodic_mode: %s', opt.aperiodic_mode);
end

out.ap.offset    = b;
out.ap.exponent  = chi;
out.ap.knee      = k;
out.ap.knee_freq = knee_freq;

out.peaks = peaks;  % [center_Hz, amplitude_log10, fwhm_Hz]

out.fit_log10.ap    = ap_fit;
out.fit_log10.peaks = peak_fit;
out.fit_log10.full  = full_fit;

out.flattened = flat_out;

out.meta = opt;
out.meta.in_range_mask = in_rng;
out.meta.noise_sd = 1.4826 * mad(flat_out(isfinite(flat_out)), 1);
out.meta.converged = converged;

% ----- full-length reconstructions on original grid -----
out.full.freq = freq(:);

% aperiodic on full grid (log10)
ff = max(out.full.freq, eps);
switch lower(opt.aperiodic_mode)
    case 'fixed'
        ap_full = b - chi * log10(ff);
    case 'knee'
        ap_full = b - log10(k + ff.^chi);
end

% peaks on full grid (log10)
peaks_full = zeros(size(ff));
for kix = 1:size(peaks,1)
    A = peaks(kix,2); mu = peaks(kix,1); fwhm = peaks(kix,3);
    sg = fwhm / (2*sqrt(2*log(2)));
    peaks_full = peaks_full + A * exp(-0.5*((ff - mu)./sg).^2);
end

out.full.ap_log10    = ap_full;
out.full.peaks_log10 = peaks_full;
out.full.model_log10 = ap_full + peaks_full;

end

% ================= helper: robust aperiodic fit with mask =================
function [params, ap_fit] = fit_aperiodic_masked(f, x, y, mode, max_clip_iters, clip_sigma, excludeMask)
% Returns:
%   fixed: params = [b, chi]
%   knee : params = [b, log10k, chi]
%   ap_fit: predicted aperiodic in log10(power) at ALL f points

mode = lower(mode);

mask = isfinite(y) & ~excludeMask;

% initial guess from linear fit in log space (fixed)
co = polyfit(x(mask), y(mask), 1);
m0 = co(1); b0 = co(2);
chi0 = max(-m0, 0.1);

% conservative knee init
log10k0 = log10(max(1e-6, (min(f(mask))+eps).^chi0 / 10));

switch mode
    case 'fixed'
        p0 = [b0, chi0];
    case 'knee'
        p0 = [b0, log10k0, chi0];
    otherwise
        error('Unknown aperiodic_mode: %s', mode);
end

% iterative sigma clipping on residuals
p = p0;
for it = 1:max_clip_iters
    yhat_all = aperiodic_eval(p, f, x, mode);

    r = y - yhat_all;
    rfit = r(mask);

    s = 1.4826 * mad(rfit, 1);
    if ~isfinite(s) || s <= 0
        s = std(rfit);
        if ~isfinite(s) || s <= 0, break; end
    end

    newmask = mask & (abs(r) <= clip_sigma*s);
    if isequal(newmask, mask)
        break;
    end
    mask = newmask;

    obj = @(pp) ap_obj(pp, f(mask), x(mask), y(mask), mode);
    p = fminsearch(obj, p, optimset('Display','off'));
end

% final refine on current mask
obj = @(pp) ap_obj(pp, f(mask), x(mask), y(mask), mode);
p = fminsearch(obj, p, optimset('Display','off'));

params = p;
ap_fit = aperiodic_eval(p, f, x, mode);

end

function yhat = aperiodic_eval(p, f, x, mode)
switch mode
    case 'fixed'
        b = p(1); chi = p(2);
        yhat = b - chi*x;
    case 'knee'
        b = p(1); log10k = p(2); chi = p(3);
        k = 10.^log10k;
        yhat = b - log10(k + f.^chi);
end
end

function J = ap_obj(p, f, x, y, mode)
% objective with soft bounds for stability
mode = lower(mode);
pen = 0;

switch mode
    case 'fixed'
        b = p(1); chi = p(2);
        yhat = b - chi*x;
        if ~isfinite(chi) || chi < 0 || chi > 8
            pen = pen + 1e6*(~isfinite(chi) + max(0,-chi) + max(0,chi-8)).^2;
        end
    case 'knee'
        b = p(1); log10k = p(2); chi = p(3);
        k = 10.^log10k;
        yhat = b - log10(k + f.^chi);

        if ~isfinite(chi) || chi < 0.1 || chi > 8
            pen = pen + 1e6*(~isfinite(chi) + max(0,0.1-chi) + max(0,chi-8)).^2;
        end
        if ~isfinite(log10k) || log10k < -12 || log10k > 12
            pen = pen + 1e6*(~isfinite(log10k) + max(0,-12-log10k) + max(0,log10k-12)).^2;
        end
end

e = y - yhat;
J = sum(e.^2) + pen;
end

% ================= helper: local maxima indices =================
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

% ================= helper: Gaussian fit cost with soft bounds =================
function J = peak_cost(p, f, r, mu_hint, sig_min, sig_max)
A  = p(1); mu = p(2); sg = p(3);
if ~isfinite(A) || ~isfinite(mu) || ~isfinite(sg) || sg <= 0
    J = 1e12; return
end
g  = A * exp(-0.5 * ((f - mu)./sg).^2);
e  = r - g;
J  = sum(e.^2);

% soft constraints
if A < 0
    J = J + 1e6 * abs(A);
end
if sg < sig_min
    J = J + 1e6 * (sig_min - sg)^2;
elseif sg > sig_max
    J = J + 1e6 * (sg - sig_max)^2;
end
% keep mu near the candidate region
if abs(mu - mu_hint) > 2*sig_max
    J = J + 1e4 * (abs(mu - mu_hint) - 2*sig_max)^2;
end
end
