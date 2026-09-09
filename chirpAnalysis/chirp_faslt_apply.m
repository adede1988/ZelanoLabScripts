function wt = chirp_faslt_apply(sig, bank)
% CHIRP_FASLT_APPLY  Apply a cached FASLT wavelet bank to one signal -> POWER [Nf x N].
%   wt = chirp_faslt_apply(sig, bank)   % sig 1xN real (or row); bank from chirp_faslt_bank
%
%   Single-buffer equivalent of Superlets/matlab-pure/nfaslt.m (geometric-mean pooling over
%   fractional superlet orders). Reuses bank.wavelets so the kernel set is built only once.
%   Output is POWER (magnitude^2 pooled), NOT analytic -- never derive phase/envelope from it.

    sig = sig(:)';
    F = bank.F; order_frac = bank.order_frac; wavelets = bank.wavelets;
    padding = bank.padding;
    Npoints = numel(sig);

    buffer = zeros(Npoints + 2*padding, 1);
    bufbegin = padding + 1; bufend = padding + Npoints;
    buffer(bufbegin:bufend) = sig;

    wt = zeros(numel(F), Npoints);
    for i_freq = 1:numel(F)
        if F(i_freq) == 0 || isempty(wavelets{i_freq,1})
            wt(i_freq,:) = 2*mean(abs(real(sig)));   % nfaslt's F==0 handling
            continue;
        end
        temp = ones(1, Npoints);
        n_wavelets = floor(order_frac(i_freq));
        for i_ord = 1:n_wavelets
            tempcx = conv(buffer, wavelets{i_freq,i_ord}, 'same');
            temp = temp .* (2 .* abs(tempcx(bufbegin:bufend)).^2)';
        end
        % fractional remainder order
        oi = ceil(order_frac(i_freq));
        if (fix(order_frac(i_freq)) ~= order_frac(i_freq)) && ~isempty(wavelets{i_freq,oi})
            exponent = order_frac(i_freq) - fix(order_frac(i_freq));
            tempcx = conv(buffer, wavelets{i_freq,oi}, 'same');
            temp = temp .* ((2 .* abs(tempcx(bufbegin:bufend)).^2)') .^ exponent;
        end
        root = 1 / order_frac(i_freq);
        wt(i_freq,:) = temp .^ root;
    end
end
