function [ITPC, frexUsed, meta] = helper_itpc_fromERP(outERP_peak, useVec, frex, fs, opts)
%HELPER_ITPC_FROMERP  Fast whole-window ITPC (phase consistency) at target freqs.
%
%   [ITPC, frexUsed, meta] = helper_itpc_fromERP(outERP_peak, useVec, frex, fs)
%   [ITPC, frexUsed, meta] = helper_itpc_fromERP(outERP_peak, useVec, frex, fs, opts)
%
% INPUTS
%   outERP_peak : [nChan x nTrial x nTime] ERPs
%   useVec      : [nTrial x 1] logical/0-1, use trials where useVec==1
%   frex        : [1 x nF] target frequencies (Hz)
%   fs          : sampling rate (Hz)
%   opts (optional struct):
%       .padFactor (default 4)     : zero-pad factor (>=1)
%       .win       (default 'hann'): taper: 'hann','hamming','blackman','none'
%       .doDetrend (default false) : detrend each ERP (time) before FFT
%       .fillNaN   (default true)  : linearly fill NaNs along time
%
% OUTPUTS
%   ITPC     : [nChan x numel(frexUsed)] ITPC magnitude at each frequency
%   frexUsed : frequencies actually used (in [0, fs/2])
%   meta     : struct with fields: nFFT, fAxis, bins, nTrialsUsed
%
% NOTES
%   Computes ITPC from the complex Fourier coefficient phase over the full window:
%     ITPC(ch,f) = | mean_trials( X / |X| ) |
%   This avoids bandpass/Hilbert edge sensitivity for short windows (e.g., 4 s).
%
% Adam Dede / helper drop-in

% ---------------- defaults ----------------
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts,'padFactor') || isempty(opts.padFactor), opts.padFactor = 4; end
if ~isfield(opts,'win')       || isempty(opts.win),       opts.win = 'hann'; end
if ~isfield(opts,'doDetrend') || isempty(opts.doDetrend), opts.doDetrend = false; end
if ~isfield(opts,'fillNaN')   || isempty(opts.fillNaN),   opts.fillNaN = true; end

% ---------------- checks ----------------
assert(ndims(outERP_peak)==3, 'outERP_peak must be [nChan x nTrial x nTime].');
[nChan, nTrial, nTime] = size(outERP_peak);
useVec = logical(useVec(:));
assert(numel(useVec)==nTrial, 'useVec length must match size(outERP_peak,2).');
frex = double(frex(:))';
fs   = double(fs);

% keep freqs in valid range
mF = (frex >= 0) & (frex <= fs/2) & isfinite(frex);
frexUsed = frex(mF);
if isempty(frexUsed)
    ITPC = nan(nChan,0);
    meta = struct('nFFT',[], 'fAxis',[], 'bins',[], 'nTrialsUsed',nnz(useVec));
    return
end

% ---------------- select + preprocess ----------------
dat = double(outERP_peak(:, useVec, :));    % [nChan x nUse x nTime]
[nChan, nUse, ~] = size(dat);

% mean-center per channel x trial (across time)
mu  = mean(dat, 3, 'omitnan');             % [nChan x nUse]
dat = dat - reshape(mu, [nChan nUse 1]);

% fill NaNs along time if requested
if opts.fillNaN && any(~isfinite(dat(:)))
    dat = fillmissing(dat, 'linear', 3, 'EndValues', 'nearest');
end

% detrend along time if requested (vectorized)
if opts.doDetrend
    tmp = reshape(permute(dat,[3 1 2]), nTime, []);  % [nTime x (nChan*nUse)]
    tmp = detrend(tmp);
    dat = permute(reshape(tmp, nTime, nChan, nUse), [2 3 1]);
end

% taper
w = local_make_window(opts.win, nTime);
dat = dat .* reshape(w, [1 1 nTime]);

% ---------------- FFT + ITPC ----------------
nFFT = 2^nextpow2(nTime * max(1, opts.padFactor));
F = fft(dat, nFFT, 3);                     % [nChan x nUse x nFFT]

% map freqs -> nearest FFT bin (single-sided)
bins = round(frexUsed * nFFT / fs) + 1;
bins = max(1, min(bins, floor(nFFT/2)+1));

Xf = F(:,:,bins);                          % [nChan x nUse x nF]
den = abs(Xf);
phaseVec = Xf ./ max(den, eps);            % unit complex vectors
ITPC = squeeze(abs(mean(phaseVec, 2, 'omitnan')));  % [nChan x nF]

% ---------------- meta ----------------
fAxis = (0:floor(nFFT/2)) * (fs/nFFT);
meta = struct();
meta.nFFT        = nFFT;
meta.fAxis       = fAxis;
meta.bins        = bins;
meta.nTrialsUsed = nUse;

end % function


% ===== local helper: window maker =====
function w = local_make_window(winName, nTime)
winName = lower(string(winName));
switch winName
    case "hann"
        w = hann(nTime);
    case "hamming"
        w = hamming(nTime);
    case "blackman"
        w = blackman(nTime);
    case "none"
        w = ones(nTime,1);
    otherwise
        error('Unknown opts.win "%s". Use hann/hamming/blackman/none.', winName);
end
w = w(:)';  % row
end