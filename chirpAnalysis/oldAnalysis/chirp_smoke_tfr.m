% chirp_smoke_tfr -- validate cached FASLT vs nfaslt, then TFR + ridge on synthetic downchirps.
setup_chirpAnalysis_paths(false);
C = chirp_config(); fs = C.fs;

% (1) parity: chirp_faslt_apply (cached) vs nfaslt (reference)
x = randn(1, 2000);
bank = chirp_faslt_bank(fs, C);
wt1 = chirp_faslt_apply(x, bank);
wt2 = nfaslt(x, fs, C.faslt.range, C.faslt.Nf, C.faslt.c1, C.faslt.ord, C.faslt.mult);
fprintf('FASLT parity: max abs rel diff = %.2e (should be ~0)\n', ...
    max(abs(wt1(:)-wt2(:))) / max(abs(wt2(:))));
fprintf('bank max half-support = %.3f s (padSec=%.2f)\n', bank.maxHalfSupportSec, C.padSec);

% (2) synthetic single downchirps 45->33 Hz, chirp at +0.2..+1.2 s post finalOnset
rng(1); nTrial = 8; blkN = 6*fs; preMs = 2;  % onset 2 s into each 6 s block
cont = []; onsets = [];
for k = 1:nTrial
    t = (0:blkN-1)/fs; t0 = preMs+0.2; t1 = preMs+1.2; f0=45; f1=33;
    inst = f0 + (f1-f0).*min(1,max(0,(t-t0)/(t1-t0)));
    ph = 2*pi*cumsum(inst)/fs;
    env = double(t>=t0 & t<=t1); env = smoothdata(env,'gaussian',round(0.08*fs));
    xb = env.*sin(ph) + 0.7*randn(1,blkN);
    onsets(end+1) = numel(cont) + preMs*fs;  %#ok<SAGROW>
    cont = [cont, xb];                         %#ok<AGROW>
end
E = chirp_epoch(cont, fs, onsets, C);
fprintf('epoch: valid=%d/%d, core=%d samp, tMs=[%.0f..%.0f]\n', ...
    sum(E.valid), nTrial, sum(E.coreIdx), E.tMs(1), E.tMs(end));
tfr = chirp_tfr_faslt(E, C);
fprintf('tfr.power [%d x %d x %d], padOK=%d\n', size(tfr.power,1),size(tfr.power,2),size(tfr.power,3),tfr.padOK);

C2 = C; C2.surr.n = 40;   % cheap surrogates for the smoke
tic; R = chirp_ridge(tfr, E, C2); el = toc;
A = R.anchors;
fprintf('ridge done in %.1fs. anchors f_hi=%.1f f_lo=%.1f (expect ~44 / ~34)\n', el, A.f_hi, A.f_lo);
for i = 1:min(4,numel(R.trial))
    T = R.trial(i); tw = T.transitionWin; if isempty(tw), tw=[NaN NaN]; end
    cw = T.coexistWin; if isempty(cw), cw=[NaN NaN]; end
    % fhat at +700 ms should be mid-chirp ~39 Hz
    f700 = NaN; if ~isempty(T.fhat), f700 = interp1(R.tMs, T.fhat, 700, 'linear', NaN); end
    fprintf('  trial %d: incl=%d nRidgesSig=%d transWin=[%.0f %.0f] coexist=[%.0f %.0f] fhat@700ms=%.1f dip=%.2f ovlp=%.2f\n', ...
        i, T.included, T.nRidgesSig, tw(1), tw(2), cw(1), cw(2), f700, T.powerDipStat, T.ridgeOverlap);
end
fprintf('SMOKE_TFR_DONE\n');
