% lib_scores.m  — shared helper functions for behavioral + macBP gamma scoring.
% (collected as local functions callable via function handles returned by lib_scores())
function H = lib_scores()
    H.welch_psd   = @welch_psd;
    H.macbp_gamma = @macbp_gamma;
    H.cue_score   = @cue_score;
    H.thresh_score= @thresh_score;
    H.o15_score   = @o15_score;
    H.spike_frac  = @spike_frac;
    H.coerce_num  = @coerce_num;
end

function [pxx,f] = welch_psd(x, fs, winSec, ovFrac, nfft)
    x = x(:); x = fillmissing(x,'linear'); x = fillmissing(x,'nearest');
    x = detrend(x,0);
    w = round(winSec*fs); if mod(w,2)==1, w=w+1; end
    win = hann(w,'periodic'); nov = round(w*ovFrac);
    [pxx,f] = pwelch(x, win, nov, nfft, fs);
end

function R = macbp_gamma(x, fs, band)
% one macBP channel -> struct with peak flattened power (dB over aperiodic) in band,
% fooof gamma-peak detectability, aperiodic exponent, spike fraction.
    R = struct('peakFlatDb',NaN,'gammaDetected',false,'peakFreq',NaN,'peakHeight',NaN, ...
               'apExponent',NaN,'flatAtPeakFreq',NaN,'spikeFrac',NaN);
    [pxx,f] = welch_psd(x, fs, 2, 0.5, 2048);
    try
        out = fooof_basic(f, pxx, 'f_range',[2 58], 'aperiodic_mode','knee', ...
                          'max_peaks',6, 'peak_width_limits',[1 12], 'peak_thresh',2);
    catch
        out = fooof_basic(f, pxx, 'f_range',[2 58], 'aperiodic_mode','fixed', ...
                          'max_peaks',6, 'peak_width_limits',[1 12], 'peak_thresh',2);
    end
    ff = f(out.meta.in_range_mask);
    flat = out.flattened;                 % log10 power over aperiodic, on fit range
    inB = ff>=band(1) & ff<=band(2);
    if any(inB)
        R.peakFlatDb = 10*max(flat(inB));  % max dB-over-aperiodic in band
    end
    R.apExponent = out.ap.exponent;
    pk = out.peaks;                        % [center height fwhm], height in log10
    if ~isempty(pk)
        ib = pk(:,1)>=band(1) & pk(:,1)<=band(2);
        if any(ib)
            R.gammaDetected = true;
            sub = pk(ib,:);
            [~,mi] = max(sub(:,2));
            R.peakFreq   = sub(mi,1);
            R.peakHeight = 10*sub(mi,2);   % dB
        end
    end
    R.spikeFrac = spike_frac(x, fs, 10, 10);   % winMs=10, K=10 sample-level
end

function sf = spike_frac(x, fs, winMs, K)
    x = x(:); x = fillmissing(x,'linear'); x = fillmissing(x,'nearest');
    w = max(2, round(winMs/1000*fs));
    d = movmax(x,w) - movmin(x,w);
    med = median(d,'omitnan');
    sigma = 1.4826*median(abs(d-med),'omitnan'); if ~isfinite(sigma)||sigma<=0, sigma=eps; end
    zd = (d-med)./sigma;
    sf = mean(zd>K);
end

function v = coerce_num(c)
    if isnumeric(c), v = double(c); return; end
    v = double(string(c));
end

function S = cue_score(bd)
    S = struct('cue_d',NaN,'cue_HR',NaN,'cue_FA',NaN,'cue_hits',NaN,'cue_misses',NaN, ...
               'cue_fas',NaN,'cue_crs',NaN,'nTrials',NaN);
    n = coerce_num(bd.n);
    [~,ia] = unique(n,'stable'); bd = bd(ia,:);
    resp = string(bd.respString);
    cue  = coerce_num(bd.cue); odor = coerce_num(bd.odor);
    keep = (resp=="Yes" | resp=="No");
    cue=cue(keep); odor=odor(keep); resp=resp(keep);
    sig = cue==odor;
    hits=sum(sig & resp=="Yes"); misses=sum(sig & resp=="No");
    fas =sum(~sig & resp=="Yes"); crs =sum(~sig & resp=="No");
    nSig=hits+misses; nNoise=fas+crs;
    S.cue_hits=hits; S.cue_misses=misses; S.cue_fas=fas; S.cue_crs=crs; S.nTrials=numel(cue);
    if nSig>0, S.cue_HR=hits/nSig; end
    if nNoise>0, S.cue_FA=fas/nNoise; end
    if nSig>0 && nNoise>0
        HRadj=(hits+0.5)/(nSig+1); FAadj=(fas+0.5)/(nNoise+1);
        S.cue_d = norminv(HRadj)-norminv(FAadj);
    end
end

function S = thresh_score(bd)
    S = struct('thresh_none',NaN,'thresh_low',NaN,'thresh_high',NaN, ...
               'thresh_low_cal',NaN,'thresh_high_cal',NaN,'nTrials',NaN);
    n = coerce_num(bd.n);
    [~,ia] = unique(n,'stable'); bd = bd(ia,:);
    lvl = coerce_num(bd.odor); inten = coerce_num(bd.intensity);
    m = @(k) mean(inten(lvl==k),'omitnan');
    S.thresh_none=m(1); S.thresh_low=m(2); S.thresh_high=m(3);
    S.thresh_low_cal  = S.thresh_low  - S.thresh_none;
    S.thresh_high_cal = S.thresh_high - S.thresh_none;
    S.nTrials = height(bd);
end

function S = o15_score(bd)
    S = struct('O15_score',NaN,'O15_acc',NaN,'nTrials',NaN);
    n = coerce_num(bd.n);
    ok = ~isnan(n);
    bd = bd(ok,:); n = n(ok);
    [~,ia] = unique(n,'stable'); bd = bd(ia,:);
    es = coerce_num(bd.expScore);
    S.O15_score = sum(es,'omitnan');
    S.O15_acc   = S.O15_score/15;
    S.nTrials   = height(bd);
end
