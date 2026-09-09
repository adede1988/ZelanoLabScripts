function P = slt_power_cont(x, fs, F, c1, ord, mult)
% SLT_POWER_CONT  Fractional adaptive (multiplicative) superlet POWER of a
% continuous 1-D signal, via cached wavelet bank + FFT convolution.
%
% Self-contained re-implementation of the superlet transform (Moca et al.
% 2021, Nat Commun; matches the Superlets repo matlab-pure/faslt.m cxmorlet
% + geometric-mean pooling exactly), adapted to run once on a whole
% recording (bank built once, conv done in the frequency domain).
%
%   P = slt_power_cont(x, fs, F, c1, ord, mult)
%     x    : 1 x N (or N x 1) real signal
%     F    : 1 x nF frequencies of interest (Hz)
%     c1   : base number of cycles (e.g. 3)
%     ord  : [oMin oMax] fractional superresolution order interval (e.g. [3 30])
%     mult : 1 multiplicative (n_cyc = i_ord*c1), 0 additive (n_cyc = i_ord+c1)
%   P    : nF x N superlet POWER (same convention as faslt output)
%
% The wavelet bank is cached (persistent) keyed on (fs,F,c1,ord,mult) so
% repeated calls within a session/batch reuse it.

    x = x(:).'; N = numel(x);
    [wav, order_frac] = get_bank(fs, F, c1, ord, mult);
    nF = numel(F);

    % FFT length: full linear convolution of x (N) with longest wavelet (L)
    Lmax = 0;
    for i=1:nF, for k=1:numel(wav{i}), Lmax = max(Lmax, numel(wav{i}{k})); end, end
    nfft = 2^nextpow2(N + Lmax - 1);
    Xf = fft(x, nfft);

    P = zeros(nF, N);
    for i = 1:nF
        of  = order_frac(i);
        nInt = floor(of);
        temp = ones(1, N);
        for k = 1:nInt
            temp = temp .* conv_same_fft(Xf, wav{i}{k}, N, nfft);
        end
        % fractional remainder wavelet (index = ceil(of))
        if of ~= fix(of)
            kf = ceil(of);
            if kf <= numel(wav{i}) && ~isempty(wav{i}{kf})
                exponent = of - fix(of);
                temp = temp .* (conv_same_fft(Xf, wav{i}{kf}, N, nfft) .^ exponent);
            end
        end
        P(i,:) = temp .^ (1/of);
    end
end

% ---- restricted (conv 'same') convolution via precomputed X fft, returns 2*|.|^2 POWER ----
function pw = conv_same_fft(Xf, w, N, nfft)
    L = numel(w);
    Wf = fft(w(:).', nfft);
    y  = ifft(Xf .* Wf);              % full linear conv, length nfft (first N+L-1 valid)
    % conv 'same' central part (matches MATLAB conv(x,w,'same')): start = floor(L/2)+1
    s  = floor(L/2) + 1;
    yc = y(s : s + N - 1);
    pw = 2 * abs(yc).^2;
end

% ---- cached bank of complex Morlet wavelets: wav{i}{k}, k=1..ceil(order_frac(i)) ----
function [wav, order_frac] = get_bank(fs, F, c1, ord, mult)
    persistent CACHE
    key = sprintf('%.6g_%s_%.6g_%.6g_%.6g_%d', fs, mat2str(F(:).',6), c1, ord(1), ord(2), mult);
    if ~isempty(CACHE) && isKey(CACHE, key)
        s = CACHE(key); wav = s.wav; order_frac = s.of; return;
    end
    nF = numel(F);
    order_frac = linspace(ord(1), ord(2), nF);
    order_int  = ceil(order_frac);
    wav = cell(nF,1);
    for i = 1:nF
        wav{i} = cell(order_int(i),1);
        for k = 1:order_int(i)
            if mult ~= 0, nc = k*c1; else, nc = k + c1; end
            wav{i}{k} = cxmorlet(F(i), nc, fs);
        end
    end
    if isempty(CACHE), CACHE = containers.Map('KeyType','char','ValueType','any'); end
    CACHE(key) = struct('wav',{wav},'of',order_frac);
end

% ---- complex Morlet, identical to Superlets/matlab-pure/faslt.m ----
function w = cxmorlet(Fc, Nc, Fs)
    sd  = (Nc/2) * (1/Fc) / 2.5;
    wl  = 2*floor(fix(6*sd*Fs)/2) + 1;
    off = fix(wl/2);
    t   = ((1:wl) - 1 - off) / Fs;
    cnorm = 1/(sd*sqrt(2*pi));
    g   = cnorm * exp(-(t.^2)/(2*sd^2));         % gaussian
    w   = exp(2i*pi*Fc*t) .* g;                  % complex morlet (bw_cf with bw=sd)
    w   = w ./ sum(g);                           % normalize by gaussian sum (as faslt)
    w   = w(:).';
end
