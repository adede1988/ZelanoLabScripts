function S = chirp_surrogates(sig, nSurr, seed)
% CHIRP_SURROGATES  Shared FFT phase-randomized surrogate generator (F6).
%   S = chirp_surrogates(sig, nSurr, seed)  -> S [nSurr x N] real surrogates
%
%   Preserves the magnitude spectrum of sig, randomizes phases (conjugate-symmetric so the
%   inverse transform is real), giving a null with identical power spectrum but destroyed
%   phase structure. ONE generator feeds ridge-2, beat modSNR, and chirplet-atom significance
%   so their 95th-percentile thresholds are mutually comparable. Deterministic given seed.

    sig = sig(:)'; N = numel(sig);
    X = fft(sig);
    rng(seed, 'twister');
    S = zeros(nSurr, N);
    half = floor(N/2);
    for s = 1:nSurr
        ph = angle(X);
        rp = 2*pi*rand(1, half) - pi;             % random phases for positive freqs (excl DC)
        ph(2:half+1) = rp;
        if mod(N,2) == 0
            ph(half+2:end) = -fliplr(rp(1:half-1)); % conjugate symmetry (Nyquist real)
        else
            ph(half+2:end) = -fliplr(rp(1:half));
        end
        ph(1) = 0;                                 % DC phase 0
        Y = abs(X) .* exp(1i*ph);
        S(s,:) = real(ifft(Y));
    end
end
