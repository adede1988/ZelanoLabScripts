function bank = chirp_faslt_bank(fs, C)
% CHIRP_FASLT_BANK  Build (and cache) the FASLT wavelet bank for reuse across trials/surrogates.
%   bank = chirp_faslt_bank(fs, C)
%
%   Replicates the wavelet-set construction inside Superlets/matlab-pure/nfaslt.m so the
%   (expensive) bank is built ONCE and reused by chirp_faslt_apply for every trial AND every
%   phase-randomized surrogate -- the batch bottleneck the spec flags (kernel caching, 6.2).
%   A persistent cache keys on (fs, range, Nf, c1, ord, mult).
%
%   bank fields: .F (1xNf Hz) .order_frac .order_int .wavelets {Nf x maxOrd} .padding
%                .maxHalfSupportSec (max wavelet half-support, for the pad assertion, F1)

    persistent CACHE
    if isempty(CACHE), CACHE = containers.Map('KeyType','char','ValueType','any'); end

    Fi   = C.faslt.range; Nf = C.faslt.Nf; c1 = C.faslt.c1; o = C.faslt.ord; mult = C.faslt.mult;
    key  = sprintf('%g_%g_%g_%d_%g_%g_%g_%g', fs, Fi(1), Fi(2), Nf, c1, o(1), o(2), mult);
    if isKey(CACHE, key), bank = CACHE(key); return; end

    if Fi(1) > Fi(2), Fi = Fi([2 1]); end
    F = linspace(Fi(1), Fi(2), Nf);
    order_frac = linspace(o(1), o(2), Nf);     % positive-only frequency domain
    order_int  = ceil(order_frac);

    wavelets = cell(numel(F), max(order_int));
    padding = 0;
    for i_freq = 1:numel(F)
        if F(i_freq) == 0, wavelets{i_freq,1} = []; continue; end
        for i_ord = 1:order_int(i_freq)
            if mult ~= 0, n_cyc = i_ord * c1; else, n_cyc = i_ord + c1; end
            wavelets{i_freq,i_ord} = cxmorlet(-F(i_freq), n_cyc, fs);   % NOTE: -F (analytic dir), as nfaslt
            padding = max(padding, fix(numel(wavelets{i_freq,i_ord})/2));
        end
    end

    bank = struct('F', F, 'order_frac', order_frac, 'order_int', order_int, ...
        'wavelets', {wavelets}, 'padding', padding, 'fs', fs, ...
        'maxHalfSupportSec', padding/fs);
    CACHE(key) = bank;
end

% --- cxmorlet + helpers: copied verbatim from Superlets/matlab-pure/nfaslt.m for parity ---
function w = cxmorlet(Fc, Nc, Fs)
    if Fc == 0, w = []; return; end
    sd  = (Nc/2) * abs(1/Fc) / 2.5;
    wl  = 2 * floor(fix(6*sd*Fs)/2) + 1;
    w   = zeros(wl,1); gi = 0; off = fix(wl/2);
    for i = 1:wl
        t    = (i-1-off)/Fs;
        w(i) = bw_cf(t, sd, Fc);
        gi   = gi + gauss(t, sd);
    end
    w = w ./ gi;
end

function res = bw_cf(t, bw, cf)
    cnorm = 1/(bw*sqrt(2*pi));
    exp1  = cnorm*exp(-(t^2)/(2*bw^2));
    res   = exp(2i*pi*cf*t)*exp1;
end

function res = gauss(t, sd)
    cnorm = 1/(sd*sqrt(2*pi));
    res   = cnorm*exp(-(t^2)/(2*sd^2));
end
