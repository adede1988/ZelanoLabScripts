function [ITPC_raw, ITPC_z, frexUsed, nullMean, nullStd, meta] = helper_itpc_timeZ_circShift(outERP, useVec, frex, fs, opts)
%HELPER_ITPC_TIMEZ_CIRCSHIFT  Time-resolved ITPC (ch x frex x time) + z via circular-shift null.
%
%   [ITPC_raw, ITPC_z, frexUsed, nullMean, nullStd, meta] = ...
%       helper_itpc_timeZ_circShift(outERP, useVec, frex, fs)
%   ... = helper_itpc_timeZ_circShift(..., opts)
%
% INPUTS
%   outERP : [nChan x nTrial x nTime]  (e.g., 32 x breaths x time)
%   useVec : [nTrial x 1] logical/0-1, use trials where useVec==1
%   frex   : [1 x nF] target freqs (Hz)
%   fs     : sampling rate (Hz)
%   opts (optional struct):
%     .nShuf         (default 100)    number of circular-shift shuffles
%     .padSec        (default 1.0)    mirror-padding seconds on each side (for wavelet convolution)
%     .nCycles       (default [])     scalar or [1 x nF] cycles for Morlet; default ramps 3..10
%     .waveletWidth  (default 3.5)    Gaussian width in SD units (time spans +/- width*s)
%     .fillNaN       (default true)   fill NaNs along time within each trial
%     .doDetrend     (default false)  detrend each trial (time) before phase extraction
%     .usePow2FFT    (default true)   use nextpow2 FFT length for convolution
%     .shuffleMethod (default 'phaseRamp')  'phaseRamp' (fast) or 'circshift' (simple)
%     .seed          (default [])     rng seed for reproducibility
%
% OUTPUTS
%   ITPC_raw : [nChan x nFused x nTime] observed ITPC(t)
%   ITPC_z   : [nChan x nFused x nTime] z-scored vs shuffle null at each timepoint
%   frexUsed : freqs used (clipped to [0, fs/2])
%   nullMean : [nChan x nFused x nTime] mean ITPC(t) across shuffles
%   nullStd  : [nChan x nFused x nTime] std  ITPC(t) across shuffles
%   meta     : struct with details (padN, nConv, halfWave, nTrialsUsed, etc.)
%
% METHOD
%   1) Extract analytic signal with Morlet wavelets (mirror padded to reduce edge artifacts)
%   2) Convert to unit phasors: U = z ./ |z|
%   3) Time-resolved ITPC: ITPC(t)=|mean_trials(U(:,t))|
%   4) Null: independently circularly shift each trial's PHASE time series by random lags,
%      recompute ITPC(t) per shuffle, then z-score per timepoint.

% ---------------- defaults ----------------
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts,'nShuf')         || isempty(opts.nShuf),         opts.nShuf = 100; end
if ~isfield(opts,'padSec')        || isempty(opts.padSec),        opts.padSec = 1.0; end
if ~isfield(opts,'nCycles'),                                   opts.nCycles = []; end
if ~isfield(opts,'waveletWidth')  || isempty(opts.waveletWidth),  opts.waveletWidth = 3.5; end
if ~isfield(opts,'fillNaN')       || isempty(opts.fillNaN),       opts.fillNaN = true; end
if ~isfield(opts,'doDetrend')     || isempty(opts.doDetrend),     opts.doDetrend = false; end
if ~isfield(opts,'usePow2FFT')    || isempty(opts.usePow2FFT),    opts.usePow2FFT = true; end
if ~isfield(opts,'shuffleMethod') || isempty(opts.shuffleMethod), opts.shuffleMethod = 'phaseRamp'; end
if ~isfield(opts,'seed'), opts.seed = []; end

assert(ndims(outERP)==3, 'outERP must be [nChan x nTrial x nTime].');
[nChan, nTrial, nTime] = size(outERP);

useVec = logical(useVec(:));
assert(numel(useVec)==nTrial, 'useVec length must match size(outERP,2).');

frex = double(frex(:))';
fs   = double(fs);

% keep freqs in [0, fs/2]
mF = isfinite(frex) & frex >= 0 & frex <= fs/2 & frex > 0;
frexUsed = frex(mF);
nF = numel(frexUsed);
if nF==0
    ITPC_raw = nan(nChan,0,nTime);
    ITPC_z   = nan(nChan,0,nTime);
    nullMean = nan(nChan,0,nTime);
    nullStd  = nan(nChan,0,nTime);
    meta = struct('nTrialsUsed',nnz(useVec));
    return
end

if ~isempty(opts.seed), rng(opts.seed); end

% cycles: scalar or vector; default ramp 3..10 cycles across freqs
if isempty(opts.nCycles)
    nCycles = linspace(3, 10, nF);
else
    nCycles = double(opts.nCycles(:))';
    if isscalar(nCycles)
        nCycles = repmat(nCycles, 1, nF);
    else
        assert(numel(nCycles)==nF, 'opts.nCycles must be scalar or length == numel(frexUsed).');
    end
end

% padding (mirror)
padN = round(opts.padSec * fs);
padN = min(padN, nTime-1);                 % safety
nTimePad = nTime + 2*padN;

% ---------------- precompute wavelet FFTs (once) ----------------
waveFFT  = cell(1,nF);
halfWave = zeros(1,nF);
nWave    = zeros(1,nF);

for fi = 1:nF
    fi
    f = frexUsed(fi);
    s = nCycles(fi) / (2*pi*f);            % Gaussian SD in seconds
    t = (-opts.waveletWidth*s : 1/fs : opts.waveletWidth*s);
    w = exp(2*1i*pi*f.*t) .* exp(-(t.^2) ./ (2*s^2));

    % normalize energy (optional but nice)
    w = w ./ sqrt(sum(abs(w).^2));

    nWave(fi) = numel(w);
    halfWave(fi) = floor(nWave(fi)/2);

    waveFFT{fi} = w; %#ok<AGROW> (temp store, FFT length set after nConv decided)
end

nConvBase = nTimePad + max(nWave) - 1;
if opts.usePow2FFT
    nConv = 2^nextpow2(nConvBase);
else
    nConv = nConvBase;
end

% finalize wavelet FFTs at nConv
for fi = 1:nF
    waveFFT{fi} = fft(waveFFT{fi}, nConv);
end

% ---------------- allocate outputs ----------------
ITPC_raw = nan(nChan, nF, nTime, 'single');
ITPC_z   = nan(nChan, nF, nTime, 'single');
nullMean = nan(nChan, nF, nTime, 'single');
nullStd  = nan(nChan, nF, nTime, 'single');

% selection
dat = double(outERP(:, useVec, :));                 % [nChan x nUse x nTime]
[nChan, nUse, ~] = size(dat);

% constants for shuffling
k = 0:(nTime-1);                                    % FFT bin indices for phase ramps
twopi_overN = 2*pi / nTime;

usePhaseRamp = strcmpi(string(opts.shuffleMethod), "phaseRamp");

% ---------------- main loops ----------------
parfor ch = 1:nChan
ch
    X = squeeze(dat(ch,:,:));                       % [nUse x nTime]
    if nUse==1, X = reshape(X, 1, nTime); end

    % mean-center per trial
    X = X - mean(X, 2, 'omitnan');

    % fill NaNs if requested
    if opts.fillNaN && any(~isfinite(X(:)))
        X = fillmissing(X, 'linear', 2, 'EndValues', 'nearest');
    end

    % detrend if requested
    if opts.doDetrend
        X = detrend(X')';                           % detrend along time for each row
    end

    % mirror padding (avoid duplicating endpoint sample)
    if padN > 0
        left  = fliplr(X(:, 2:padN+1));
        right = fliplr(X(:, end-padN:end-1));
        Xpad  = [left, X, right];                   % [nUse x nTimePad]
    else
        Xpad  = X;
    end

    Xfft = fft(Xpad, nConv, 2);                     % [nUse x nConv]

    for fi = 1:nF
        % convolution -> analytic signal
        convRes = ifft(Xfft .* waveFFT{fi}, [], 2); % [nUse x nConv]
        hw = halfWave(fi);

        % crop to padded length then to original epoch
        convRes = convRes(:, hw+1 : hw+nTimePad);   % [nUse x nTimePad]
        if padN > 0
            convRes = convRes(:, padN+1 : padN+nTime); % [nUse x nTime]
        end

        % unit phasors (phase time series encoded as complex unit vectors)
        U = convRes ./ max(abs(convRes), eps);      % [nUse x nTime], complex

        % observed ITPC(t)
        itpcObs = abs(mean(U, 1));                  % [1 x nTime]
        ITPC_raw(ch, fi, :) = single(itpcObs);

        % ---------- shuffle null ----------
        nShuf = max(1, round(opts.nShuf));
        muNull = zeros(1, nTime, 'double');
        M2     = zeros(1, nTime, 'double');

        if usePhaseRamp
            % fast: shift in frequency domain of unit phasors
            Uf = fft(U, [], 2);                     % [nUse x nTime]
            for ss = 1:nShuf
                sh = randi([0 nTime-1], [nUse 1]);  % random lag per trial

                % phase ramp per trial (vectorized). ramp size [nUse x nTime]
                ramp = exp(-1i * twopi_overN * (double(sh) * k));
                meanF = mean(Uf .* ramp, 1);        % [1 x nTime]
                meanU = ifft(meanF, [], 2);         % [1 x nTime]
                itpcSh = abs(meanU);                % [1 x nTime]

                if ss == 1
                    muNull = itpcSh;
                else
                    d = itpcSh - muNull;
                    muNull = muNull + d / ss;
                    M2 = M2 + d .* (itpcSh - muNull);
                end
            end

        else
            % simple: shift in time domain (no big exponentials, but loops over trials)
            for ss = 1:nShuf
                sh = randi([0 nTime-1], [nUse 1]);
                sumShift = zeros(1, nTime, 'like', U(1,:));
                for tr = 1:nUse
                    sumShift = sumShift + circshift(U(tr,:), sh(tr), 2);
                end
                itpcSh = abs(sumShift ./ nUse);

                itpcSh = double(itpcSh);
                if ss == 1
                    muNull = itpcSh;
                else
                    d = itpcSh - muNull;
                    muNull = muNull + d / ss;
                    M2 = M2 + d .* (itpcSh - muNull);
                end
            end
        end

        if nShuf > 1
            sdNull = sqrt(M2 / (nShuf - 1));
        else
            sdNull = nan(size(muNull));
        end

        nullMean(ch, fi, :) = single(muNull);
        nullStd(ch, fi, :)  = single(sdNull);

        % z-score per timepoint (guard tiny std)
        sdNull(sdNull < eps) = NaN;
        itpcZ = (double(itpcObs) - muNull) ./ sdNull;

        ITPC_z(ch, fi, :) = single(itpcZ);
    end
end

% ---------------- meta ----------------
meta = struct();
meta.fs          = fs;
meta.nTrialsUsed = nUse;
meta.nTime       = nTime;
meta.padN        = padN;
meta.nTimePad    = nTimePad;
meta.nConv       = nConv;
meta.frexUsed    = frexUsed;
meta.nCycles     = nCycles;
meta.halfWave    = halfWave;
meta.nShuf       = max(1, round(opts.nShuf));
meta.shuffleMethod = char(string(opts.shuffleMethod));

end