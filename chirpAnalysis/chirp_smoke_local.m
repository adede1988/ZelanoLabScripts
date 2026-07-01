% _smoke_local -- confirm the chirp engine chain runs (home or lab). Synthetic only.
%   matlab -batch "run('...\chirpAnalysis\_smoke_local.m')"
setup_chirpAnalysis_paths(true);
C = chirp_config(); fs = C.fs;
fprintf('\n--- engines on path: nfaslt=%d tfridge=%d mp_adapt_chirplets=%d hilbert=%d\n', ...
    ~isempty(which('nfaslt')), ~isempty(which('tfridge')), ...
    ~isempty(which('mp_adapt_chirplets')), ~isempty(which('hilbert')));

% synthetic single downchirp 46->32 Hz over [1,2] s within a 4 s trial, + noise
t = (0:fs*4-1)/fs; tc0=1; tc1=2; f0=46; f1=32;
inst = f0 + (f1-f0).*min(1,max(0,(t-tc0)/(tc1-tc0)));
ph = 2*pi*cumsum(inst)/fs;
env = double(t>=tc0 & t<=tc1); env = smoothdata(env,'gaussian',round(0.1*fs));
x = env.*sin(ph) + 0.5*randn(size(t));

% --- FASLT ---
try
    wt = nfaslt(x, fs, C.faslt.range, C.faslt.Nf, C.faslt.c1, C.faslt.ord, C.faslt.mult);
    F  = linspace(C.faslt.range(1), C.faslt.range(2), C.faslt.Nf);
    fprintf('FASLT ok: size [%d x %d], F=[%.0f..%.0f]\n', size(wt,1), size(wt,2), F(1), F(end));
catch ME, fprintf('FASLT FAIL: %s\n', ME.message); end

% --- tfridge orientation test (rows=freq) ---
try
    inB = F>=C.ridge.band(1) & F<=C.ridge.band(2);
    sst = wt(inB,:); fb = F(inB);
    [fr,ir,lr] = tfridge(sst, fb, C.ridge.penalty);   %#ok<ASGLU>
    fprintf('tfridge ok (rows=freq): fr len=%d, fr(0.5s)=%.1f fr(1.5s)=%.1f Hz\n', ...
        numel(fr), fr(round(0.5*fs)), fr(round(1.5*fs)));
catch ME
    fprintf('tfridge rows=freq FAIL: %s -- trying transpose\n', ME.message);
    try, [fr,ir,lr] = tfridge(sst.', fb, C.ridge.penalty); fprintf('tfridge ok (transpose)\n');
    catch ME2, fprintf('tfridge transpose FAIL: %s\n', ME2.message); end
end

% --- MPACT chirplet on analytic of bandpassed burst over the transition window ---
try
    bb = filtfilt(fir1(round(0.2*fs)*2, C.bbBand/(fs/2)), 1, x);  % crude bandpass for smoke
    seg = bb(t>=0.9 & t<=2.1);                                    % transition window
    xa = hilbert(seg);
    Q=2; M=C.chirplet.M; D=C.chirplet.D; i0=C.chirplet.i0; rad=C.chirplet.radix; mn=C.chirplet.mnits;
    P = mp_adapt_chirplets(double(xa(:)), Q, M, D, i0, rad, 'no', mn, 2, 'RefineAlgorithm','expectmax','PType','Oneill');
    fcHz = P(:,3).*fs/(2*pi); crHzs = P(:,4).*fs^2/(2*pi); tcS = P(:,2)/fs;
    fprintf('MPACT ok: %d atoms\n', size(P,1));
    for k=1:size(P,1)
        fprintf('  atom%d |A|=%.2f tc=%.3fs fc=%.1fHz chirp=%.1fHz/s d=%.1fsamp\n', ...
            k, abs(P(k,1)), tcS(k), fcHz(k), crHzs(k), P(k,5));
    end
catch ME, fprintf('MPACT FAIL: %s\n  %s\n', ME.message, getReport(ME,'basic')); end

% --- chirp_epoch + chirp_burstlen sanity ---
try
    cont = repmat(x, 1, 3);                 % 12 s "continuous"
    ons  = [2*fs, 6*fs, 10*fs] + tc0*fs;    % onsets so the chirp sits at +0..+1s post-onset
    E = chirp_epoch(cont, fs, ons, C);
    fprintf('chirp_epoch ok: dataPad [%d x %d], core %d samp, valid=%d/%d\n', ...
        size(E.dataPad,1), size(E.dataPad,2), sum(E.coreIdx), sum(E.valid), numel(ons));
    % gamma power per epoch (broadband env^2), core-trimmed
    gp = nan(numel(ons), sum(E.coreIdx));
    for i=1:numel(ons)
        if ~E.valid(i), continue; end
        bbp = filtfilt(fir1(round(0.2*fs)*2, C.bbBand/(fs/2)),1,E.dataPad(i,:));
        envp = abs(hilbert(bbp)).^2; gp(i,:) = envp(E.coreIdx);
    end
    B = chirp_burstlen(gp, E.tMs, C);
    fprintf('chirp_burstlen ok: thr=%.3f nBurst=%d/%d median len=%.0f ms\n', ...
        B.thr, B.summary.nBurst, B.summary.nValid, B.summary.median);
catch ME, fprintf('epoch/burst FAIL: %s\n', ME.message); end

fprintf('SMOKE_DONE\n');
