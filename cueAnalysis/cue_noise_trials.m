function out = cue_noise_trials(sig, fs, events, epWin, threshUV, winMs)
% CUE_NOISE_TRIALS  Per-trial sharp-deflection noise detection (analysis3).
%
%   out = cue_noise_trials(sig, fs, events [,epWin,threshUV,winMs])
%
%   Epochs the channel around each event (default [-1.75 +5.75] s, the
%   single-trial-plot window) and flags a trial as NOISY if anywhere in the
%   epoch the max-to-min range within a sliding ~winMs window exceeds threshUV
%   uV (default: >80 uV in a 10 ms window). This is the noise rule used for both
%   the red flagging on singleTrialRawMac.png AND trial rejection in the
%   spectrograms (it replaces the old z-score QC).
%
%   out fields: .D [pnts x N] (NaN cols = out of bounds) .ok [Nx1] .noisy [Nx1]
%     .times (ms) .epWin .threshUV .winMs .win

    if nargin < 4 || isempty(epWin),    epWin = [-1.75 5.75]; end
    if nargin < 5 || isempty(threshUV), threshUV = 80; end
    if nargin < 6 || isempty(winMs),    winMs = 10; end

    s0 = round(epWin(1)*fs); s1 = round(epWin(2)*fs); nF = s1 - s0 + 1;
    T = numel(sig); N = numel(events);
    D = nan(nF, N); ok = false(N,1);
    for i = 1:N
        e = events(i);
        if ~isfinite(e) || e <= 0, continue; end
        a = round(e) + s0; b = round(e) + s1;
        if a < 1 || b > T, continue; end
        seg = sig(a:b);
        if any(~isfinite(seg)), continue; end
        D(:,i) = seg(:); ok(i) = true;
    end

    win = max(2, round(winMs/1000*fs));          % samples spanning ~winMs (10 ms = 5 @500Hz)
    noisy = false(N,1);
    for i = find(ok)'
        x = D(:,i);
        rng = movmax(x, win) - movmin(x, win);   % local max-min over a winMs window
        noisy(i) = any(rng > threshUV);
    end

    out = struct('D', D, 'ok', ok, 'noisy', noisy, 'times', ((s0:s1)/fs*1000), ...
                 'epWin', epWin, 'threshUV', threshUV, 'winMs', winMs, 'win', win);
end
