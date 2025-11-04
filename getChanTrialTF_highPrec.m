function trialTF = getChanTrialTF_highPrec(trialDat, frex, nCyclesOrStd, srate, varargin)
% trialDat: [T x Ntr] real data (time x trials)
% frex    : [1 x F] frequencies (Hz)
% nCyclesOrStd: either scalar/vec of (#cycles) OR time-std (s). See 'mode'.
% srate   : sampling rate (Hz)
%
% Options:
%   'WidthMode'  : 'cycles' (default) or 'std'  (seconds)
%   'TimeWin'    : 2.0       % seconds total half-width for wavelet window (±TimeWin)
%   'UseParfor'  : true/false
%
% Returns:
%   trialTF: [T x Ntr x F] complex TF (analytic wavelet coefficients)

% ----- options -----
p = inputParser;
p.addParameter('WidthMode','cycles');
p.addParameter('TimeWin', 1.0);          % seconds each side (total window = 2*TimeWin)
p.addParameter('UseParfor', false);
p.parse(varargin{:});
opt = p.Results;

[T, Ntr]  = size(trialDat);
F         = numel(frex);
frex      = frex(:).';                   % row

% ----- mirror-pad (simple, known shape: 3*T x Ntr) -----
padDat = [flipud(trialDat); trialDat; flipud(trialDat)];

% ----- time vector in *seconds*, consistent with srate -----
time = -opt.TimeWin : 1/srate : opt.TimeWin;
n_wavelet            = numel(time);
n_data               = numel(padDat);
n_convolution        = n_wavelet + n_data - 1;
n_fft                = 2^nextpow2(n_convolution);
half_wave            = (n_wavelet - 1)/2;
if mod(n_wavelet,2)==0
    error('Wavelet length must be odd; pick TimeWin so that 2*TimeWin*srate is even.');
end

% ----- FFT(data) once -----
fftDat = fft(reshape(padDat, 1, []), n_fft);

% ----- frequency-dependent time SD (sigma_t) -----
switch lower(opt.WidthMode)
    case 'cycles'
        nC = nCyclesOrStd(:).';
        if isscalar(nC), nC = repmat(nC, 1, F); end
        sigma_t = nC ./ (2*pi*frex);          % seconds
    case 'std'
        sigma_t = nCyclesOrStd(:).';
        if isscalar(sigma_t), sigma_t = repmat(sigma_t, 1, F); end
    otherwise
        error('WidthMode must be ''cycles'' or ''std''.');
end

% ----- output (complex) -----
trialTF = complex(zeros(T, Ntr, F, 'like', trialDat));

% ----- per-frequency loop (parfor optional) -----
if opt.UseParfor
    parfor fi = 1:F
        trialTF(:,:,fi) = do_one_freq(fi, frex, sigma_t, time, n_fft, ...
                                      fftDat, n_convolution, half_wave, ...
                                      T, padDat);  % use constFft.Value / constPad.Value if using Constants
    end
else
    for fi = 1:F
        trialTF(:,:,fi) = do_one_freq(fi, frex, sigma_t, time, n_fft, ...
                                      fftDat, n_convolution, half_wave, ...
                                      T, padDat);
    end
end

end

function outSlice = do_one_freq(fi, frex, sigma_t, time, n_fft, ...
                                fftDat, n_convolution, half_wave, ...
                                T, padDat)
    f = frex(fi);
    s = sigma_t(fi);

    % Morlet, unit L2 energy
    gauss = exp(-(time.^2)/(2*s^2));
    cmplx = exp(1i*2*pi*f.*time);
    w     = gauss .* cmplx;
    w     = w ./ sqrt(sum(abs(w).^2));

    W        = fft(w, n_fft);
    convsig  = ifft(W .* fftDat);
    convsig  = convsig(1:n_convolution);
    convsig  = convsig(half_wave+1:end-half_wave);
    tmp      = reshape(convsig, size(padDat));     % [3T x Ntr]
    outSlice = tmp(T+1:2*T, :);                    % middle (unpadded) block
end