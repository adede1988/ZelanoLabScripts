function A = chirp_ridge_anchors(tfr, C)
% CHIRP_RIDGE_ANCHORS  Rough f_hi/f_lo from the trial-AVERAGED power TFR (spec 6.4).
%   A = chirp_ridge_anchors(tfr, C)  -> A.f_hi, A.f_lo (Hz, clamped C.ridge.anchorClamp)
%
%   Robust to BOTH regimes:
%     - DUAL (two coexisting tones): the burst-marginal spectrum of meanPower has two prominent
%       peaks -> f_hi/f_lo = the two tones (so fBeat=|f_hi-f_lo| is the real gap).
%     - SINGLE (sweep): the marginal is unimodal -> fall back to the primary-ridge early/late
%       endpoints (f_hi @20th-pctile transition time, f_lo @80th), i.e. the sweep span.
%   Deliberately rough; ANNOTATION / spatial-band / beat-fBeat / figures only -- never a
%   per-trial window gate, never the chirplet.

    f = tfr.freqs(tfr.bandMask); tMs = tfr.tMs; fs = tfr.fs;
    A = struct('f_hi', NaN, 'f_lo', NaN, 'src','none');
    M = tfr.meanPower(tfr.bandMask, :);
    if all(~isfinite(M(:))) || all(M(:)==0), return; end

    % burst-present time span from the time-marginal power
    pwrT = sum(M,1,'omitnan'); post = tMs>0 & tMs<=C.burst.searchMs(2);
    idx = find(post); if isempty(idx), return; end
    [pk,kr] = max(pwrT(idx)); thr = 0.5*pk; kpk = idx(kr);
    a=kpk; while a>1 && post(a-1) && pwrT(a-1)>thr, a=a-1; end
    b=kpk; while b<numel(pwrT) && post(b+1) && pwrT(b+1)>thr, b=b+1; end

    % (1) two-peak frequency marginal over the burst
    marg = mean(M(:,a:b),2,'omitnan');
    fHiM = NaN; fLoM = NaN;
    try
        [~,locs,~,prom] = findpeaks(marg,'SortStr','descend');
    catch, locs=[]; prom=[]; end
    if ~isempty(locs)
        keep = prom > 0.2*max(prom); locs = locs(keep);
        fp = f(locs);
        if numel(fp) >= 2
            f1 = fp(1); sep = abs(fp(2:end)-f1); j = find(sep>=3,1,'first');
            if ~isempty(j), pr = sort([f1 fp(1+j)]); fLoM=pr(1); fHiM=pr(2); end
        end
    end

    % (2) primary-ridge early/late endpoints (sweep span)
    [fhat,iR] = tfridge(M, f, C.ridge.penalty, 'NumRidges',1);
    fhat=fhat(:)'; iR=iR(:)'; rp = M(sub2ind(size(M),iR,1:numel(tMs)));
    smW = max(1,round(C.ridge.fhatSmoothMs/1000*fs)); fhatS = movmedian(fhat,smW);
    tA = tMs(a); tB = tMs(b);
    fHiR = interp1(tMs, fhatS, tA+0.2*(tB-tA), 'linear','extrap');
    fLoR = interp1(tMs, fhatS, tA+0.8*(tB-tA), 'linear','extrap');

    if isfinite(fHiM) && isfinite(fLoM)
        A.f_hi=fHiM; A.f_lo=fLoM; A.src='marginal2peak';
    else
        A.f_hi=max(fHiR,fLoR); A.f_lo=min(fHiR,fLoR); A.src='ridgeEndpoints';
    end
    cl = C.ridge.anchorClamp;
    A.f_hi = min(max(A.f_hi,cl(1)),cl(2));
    A.f_lo = min(max(A.f_lo,cl(1)),cl(2));
end
