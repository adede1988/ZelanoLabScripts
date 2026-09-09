function B = chirp_burstlen(gpow, tMs, C)
% CHIRP_BURSTLEN  Burst-length distribution from sniff-locked gamma power (D4, user spec).
%   B = chirp_burstlen(gpow, tMs, C)
%
%   gpow : [nTrial x nT] smoothed-or-raw gamma POWER, finalOnset-locked (e.g. broadband
%          25-58 Hz envelope^2, or FASLT band-mean). Rows of all-NaN are skipped.
%   tMs  : [1 x nT] time (ms) relative to finalOnset.
%   C    : chirp_config (uses C.burst.smoothMs/.peakPct/.searchMs).
%
%   Method (exactly as specified): smooth each trial's gamma-power series; take each trial's
%   PEAK within (0, searchMax] ms; the THRESHOLD is the C.burst.peakPct-th percentile of the
%   cross-trial peak distribution. A trial's burst spans from where its smoothed signal first
%   REACHES the threshold to where it FALLS back BELOW it (searching only within the post-onset
%   window). Trials whose peak never reaches the threshold have no measurable burst (NaN).
%
%   Fields: B.thr B.peak[nTrial] B.peakLatMs[nTrial] B.startMs[nTrial] B.stopMs[nTrial]
%     B.lenMs[nTrial] B.truncated[nTrial] B.hasBurst[nTrial] B.coverage
%     B.summary (median/p10/p25/p75/p90/max len among bursts, nBurst)

    if nargin < 3 || isempty(C), C = chirp_config(); end
    [N, nT] = size(gpow);

    dt = median(diff(tMs));                         % ms/sample
    w  = max(1, round(C.burst.smoothMs / dt));
    sm = movmean(gpow, w, 2, 'omitnan');            % smoothed power per trial

    inSearch = tMs > C.burst.searchMs(1) & tMs <= C.burst.searchMs(2);
    si = find(inSearch);

    peak = nan(N,1); peakLatMs = nan(N,1);
    valid = false(N,1);
    for i = 1:N
        if all(isnan(sm(i,:))), continue; end
        valid(i) = true;
        seg = sm(i, si);
        [pk, k] = max(seg);
        peak(i) = pk; peakLatMs(i) = tMs(si(k));
    end

    thr = prctile(peak(valid), C.burst.peakPct);    % cross-trial 75th pctile of per-trial peaks

    startMs = nan(N,1); stopMs = nan(N,1); lenMs = nan(N,1);
    truncated = false(N,1); hasBurst = false(N,1);
    for i = find(valid)'
        seg = sm(i, si); t = tMs(si);
        above = seg >= thr;
        if ~any(above), continue; end              % peak below threshold -> no burst
        hasBurst(i) = true;
        kStart = find(above, 1, 'first');
        % find fall-back-below AFTER the trial's own peak
        [~, kPk] = max(seg);
        kStop = find(~above(kPk:end), 1, 'first');
        if isempty(kStop)
            kStop = numel(seg); truncated(i) = true; % never fell below within the search window
        else
            kStop = kPk + kStop - 1;
        end
        startMs(i) = t(kStart); stopMs(i) = t(kStop);
        lenMs(i)  = stopMs(i) - startMs(i);
    end

    L = lenMs(hasBurst);
    summary = struct('nValid', sum(valid), 'nBurst', sum(hasBurst), ...
        'median', medsafe(L,50), 'p10', medsafe(L,10), 'p25', medsafe(L,25), ...
        'p75', medsafe(L,75), 'p90', medsafe(L,90), 'max', maxsafe(L), ...
        'fracTruncated', sum(truncated)/max(1,sum(hasBurst)));

    B = struct('thr', thr, 'peak', peak, 'peakLatMs', peakLatMs, 'startMs', startMs, ...
        'stopMs', stopMs, 'lenMs', lenMs, 'truncated', truncated, 'hasBurst', hasBurst, ...
        'valid', valid, 'coverage', sum(hasBurst)/max(1,sum(valid)), 'summary', summary, ...
        'peakPct', C.burst.peakPct, 'smoothMs', C.burst.smoothMs, 'searchMs', C.burst.searchMs);
end

function v = medsafe(x,p), if isempty(x), v = NaN; else, v = prctile(x,p); end, end
function v = maxsafe(x),   if isempty(x), v = NaN; else, v = max(x);       end, end
