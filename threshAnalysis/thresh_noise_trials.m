function out = thresh_noise_trials(sig, fs, events, epWin, K, winMs)
% THRESH_NOISE_TRIALS  Per-trial sharp-deflection noise detection (RELATIVE rule).
%
%   out = thresh_noise_trials(sig, fs, events [,epWin,K,winMs])
%
%   Detects sharp deflections RELATIVE to the channel's own variability, so that
%   a uniformly high-amplitude channel is not penalized for being high-amplitude
%   (the earlier absolute >80 uV/10 ms rule rejected such channels wholesale).
%
%   Method: d = max-min within a sliding ~winMs window (default 10 ms). Over the
%   WHOLE recording, robust-z each window: zd = (d - median(d)) / (1.4826*MAD(d)).
%   A window is an "event" if zd > K (it is K robust-SDs above this channel's
%   TYPICAL 10 ms swing). A trial is NOISY if its epoch contains any event.
%   K is calibrated once (thresh_calibrate_noise) to drop >0% and <=20% of trials
%   dataset-wide; the SAME K applies to every session.
%
%   The thresh anchor is TTL.start (= sniff - 2 s @500 Hz); the epoch window below
%   is centered so it spans the pre-start baseline through the post-sniff response.
%
%   out fields: .D [pnts x N] (NaN cols = OOB) .ok .noisy .flagMask [pnts x N]
%     (zd>K per sample, for plotting) .Z (zd per sample) .times(ms) .K .win .med .sigma

    if nargin < 4 || isempty(epWin), epWin = [-1.75 5.75]; end
    if nargin < 5 || isempty(K),     K = 10; end   % CALIBRATED (thresh_calibrate_noise);
    %   ported default from cue (K=10 -> ~8% dataset-wide loss). Re-run
    %   thresh_calibrate_noise on thresh amplitudes and set the final K here.
    if nargin < 6 || isempty(winMs), winMs = 10; end

    win = max(2, round(winMs/1000*fs));            % samples spanning ~winMs (10 ms = 5 @500Hz)
    d = movmax(sig, win) - movmin(sig, win);       % 10 ms max-min over the whole recording
    med = median(d, 'omitnan');
    sigma = 1.4826 * median(abs(d - med), 'omitnan');
    if ~(sigma > 0), sigma = eps; end
    zd = (d - med) ./ sigma;                       % robust z of each window

    s0 = round(epWin(1)*fs); s1 = round(epWin(2)*fs); nF = s1 - s0 + 1;
    T = numel(sig); N = numel(events);
    D = nan(nF, N); Z = nan(nF, N); ok = false(N,1);
    for i = 1:N
        e = events(i);
        if ~isfinite(e) || e <= 0, continue; end
        a = round(e)+s0; b = round(e)+s1;
        if a < 1 || b > T, continue; end
        seg = sig(a:b);
        if any(~isfinite(seg)), continue; end
        D(:,i) = seg(:); Z(:,i) = zd(a:b)'; ok(i) = true;
    end

    flagMask = Z > K;                              % per-sample relative-extreme flag
    noisy = false(N,1);
    for i = find(ok)', noisy(i) = any(flagMask(:,i)); end

    out = struct('D', D, 'ok', ok, 'noisy', noisy, 'flagMask', flagMask, 'Z', Z, ...
                 'times', ((s0:s1)/fs*1000), 'epWin', epWin, 'K', K, 'win', win, ...
                 'med', med, 'sigma', sigma);
end
