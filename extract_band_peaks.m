function bpTbl = extract_band_peaks(res_low, res_high, bands)
% Return one peak per band as a table.
% Columns: band, freq_Hz, amp_log10, fwhm_Hz, found, default_used

names = fieldnames(bands);
n = numel(names);

band      = strings(n,1);
freq_Hz   = nan(n,1);
amp_log10 = nan(n,1);
fwhm_Hz   = nan(n,1);
found     = false(n,1);
default_used = false(n,1);

for i = 1:n
    nm  = names{i};
    rng = bands.(nm);
    band(i) = string(nm);

    % Default (geometric mean of edges)
    defaultHz = 10.^mean(log10([max(rng(1), eps), rng(2)]));

    % Gather peaks from both fits that land in the band
    P = [];
    if ~isempty(res_low.peaks)
        sel = res_low.peaks(:,1) >= rng(1) & res_low.peaks(:,1) <= rng(2);
        P = [P; res_low.peaks(sel,:)];
    end
    if ~isempty(res_high.peaks)
        sel = res_high.peaks(:,1) >= rng(1) & res_high.peaks(:,1) <= rng(2);
        P = [P; res_high.peaks(sel,:)];
    end

    if ~isempty(P)
        % Pick the largest (by log-amplitude, col 2)
        [~, k]     = max(P(:,2));
        freq_Hz(i) = P(k,1);
        amp_log10(i) = P(k,2);
        fwhm_Hz(i)   = P(k,3);
        found(i)       = true;
        default_used(i)= false;
    else
        % No peak in band -> use default center
        freq_Hz(i)     = defaultHz;
        amp_log10(i)   = NaN;
        fwhm_Hz(i)     = NaN;
        found(i)       = false;
        default_used(i)= true;
    end
end

bpTbl = table(categorical(band), freq_Hz, amp_log10, fwhm_Hz, found, default_used, ...
    'VariableNames', {'band','freq_Hz','amp_log10','fwhm_Hz','found','default_used'});
end
