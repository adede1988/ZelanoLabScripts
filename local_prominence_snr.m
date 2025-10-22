function [promSNR, localMAD] = local_prominence_snr(sig, left, right, flankHalf, polarity, globalMAD, madFloorFrac)
% sig        : 1 x T vector (single channel & trial)
% left/right : indices bracketing the event (inclusive)
% flankHalf  : half-width (samples) for each flank window
% polarity   : 'pos' or 'neg' (what kind of event we test)
% globalMAD  : scalar robust scale for this channel (optional but recommended)
% madFloorFrac : fraction of globalMAD to use as a floor (e.g., 0.3)
%
% Outputs:
%   promSNR  : (center - flankRef) / max(localMAD, madFloorFrac*globalMAD, eps), for 'pos'
%              (flankRef - center) / ...                                     , for 'neg'
%   localMAD : robust scale estimated from flank samples only

T = numel(sig);
if nargin < 6 || isempty(globalMAD), globalMAD = []; end
if nargin < 7 || isempty(madFloorFrac), madFloorFrac = 0.3; end

% Flank indices (clip to bounds)
LL1 = max(1, left - flankHalf);   LL2 = left;
RR1 = right;                      RR2 = min(T, right + flankHalf);

% Collect flank samples (exclude center span as much as possible)
flankSamples = [];
if LL1 <= LL2, flankSamples = [flankSamples, sig(LL1:LL2)]; end
if RR1 <= RR2, flankSamples = [flankSamples, sig(RR1:RR2)]; end

% Robust local scale from flanks
if isempty(flankSamples)
    localMAD = NaN;
else
    medF = median(flankSamples);
    localMAD = 1.4826 * median(abs(flankSamples - medF));
end

% Floor the scale with a fraction of globalMAD (if provided)
scale = localMAD;
if ~isempty(globalMAD)
    scale = max(scale, madFloorFrac * globalMAD);
end
scale = max(scale, eps);

% Center and flank "reference" extremes with matching polarity
if strcmpi(polarity,'pos')
    center = max(sig(left:right));
    % worst (largest) positive in each flank
    Lmax = max(sig(LL1:LL2)); if isempty(Lmax), Lmax = -inf; end
    Rmax = max(sig(RR1:RR2)); if isempty(Rmax), Rmax = -inf; end
    flankRef = max([Lmax, Rmax]);
    delta = center - flankRef;  % >= 0 if center towers above flanks
else % 'neg'
    center = min(sig(left:right));
    Lmin = min(sig(LL1:LL2)); if isempty(Lmin), Lmin = +inf; end
    Rmin = min(sig(RR1:RR2)); if isempty(Rmin), Rmin = +inf; end
    flankRef = min([Lmin, Rmin]);
    delta = flankRef - center;  % >= 0 if center dips below flanks
end

promSNR = delta / scale;
end
