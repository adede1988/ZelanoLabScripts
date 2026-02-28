function [pacOut, meta] = pac_breathTemplate_timeResolvedPAC(rawDat, macRaw, keyBreathIDX, gamMed, fs, bpHz, PACfrex, varargin)
% pac_breathTemplate_timeResolvedPAC
%
% Time-resolved PAC sampled at breath-template indices.
%
% TWO MODES (selected by whether arg2 = macRaw is empty):
%   (1) Legacy (macRaw provided):
%       - Phase source: rawDat.data
%       - Amp source:   macRaw.data  (bandpass bpHz, Hilbert envelope)
%
%   (2) Resp mode (macRaw = []):
%       - Phase source: rawDat.rsp   (detrended)
%       - Amp source:   rawDat.data  (bandpass bpHz, Hilbert envelope)
%
% Phase extraction:
%   - For f < sosBelowHz: IIR bandpass around f (Butterworth), then Hilbert phase
%   - For f >= sosBelowHz: Morlet wavelet phase (waveCycles cycles), then resample complex series
%
% Window:
%   - nCyclesWin cycles at phase frequency, capped at winCapSec (default 50 s)
%
% Nulls:
%   (A) GLOBAL circular shift: one shift per shuffle applied to all windows (existing)
%   (B) BREATH-wise circular shift: one shift per breath per shuffle (new)
%
% Output:
%   pacOut: [nBreaths x nTpl x nFrex x 9]
%       (:,:,:,1) = coupling strength: |<A*u>| / <A>
%       (:,:,:,2) = corr-style:        |<A*u>| / sqrt(<A^2>)
%       (:,:,:,3) = z-score of (1) vs GLOBAL circular-shift null
%       (:,:,:,4) = z-score of (2) vs GLOBAL circular-shift null
%       (:,:,:,5) = preferred phase:   angle(<A*u>)
%       (:,:,:,6) = z-score of (1) vs BREATH-wise circular-shift null
%       (:,:,:,7) = z-score of (2) vs BREATH-wise circular-shift null
%       (:,:,:,8) = debiased coupling: |< (A-<A>) * u >| / sqrt(<(A-<A>)^2>)
%       (:,:,:,9) = z-score of (8) vs GLOBAL circular-shift null
%
% Name-value options:
%   'nCyclesWin'     (default 20)
%   'nShuf'          (default 200)
%   'winShape'       (default 'gaussian')  % 'gaussian' or 'boxcar'
%   'waveCycles'     (default 6)
%   'envLowpassHz'   (default 20)
%   'targetFs'       (default 100)         % resample both phase+amp streams to this
%   'minShiftSec'    (default 10)
%   'rngSeed'        (default 0)
%   'chunkSize'      (default 2000)
%
%   Hybrid cutoff + SOS band controls:
%   'sosBelowHz'     (default 2)
%   'bwFrac'         (default 0.30)        % bandwidth = max(bwMin, bwFrac*f)
%   'bwMin'          (default 0.06)        % Hz
%   'sosOrder'       (default 4)           % Butterworth order for bandpass (SOS)
%   'sosEdgeCycles'  (default 3)           % edge clearance (cycles of f) for SOS method
%   'winCapSec'      (default 50)          % hard cap on window length
%
% Notes:
%   - Uses resample() (anti-aliasing) rather than integer decimation.
%   - Requires Signal Processing Toolbox (resample, butter, hilbert, filtfilt).

% -------------------- parse options --------------------
p = inputParser;
p.addParameter('nCyclesWin',    20, @(x)isscalar(x)&&x>0);
p.addParameter('nShuf',        200, @(x)isscalar(x)&&x>=20);
p.addParameter('winShape', 'gaussian', @(s)ischar(s)||isstring(s));
p.addParameter('waveCycles',     6, @(x)isscalar(x)&&x>0);
p.addParameter('envLowpassHz',  20, @(x)isscalar(x)&&x>0);
p.addParameter('targetFs',     100, @(x)isscalar(x)&&x>0);
p.addParameter('minShiftSec',   10, @(x)isscalar(x)&&x>0);
p.addParameter('rngSeed',        0, @(x)isscalar(x));
p.addParameter('chunkSize',   2000, @(x)isscalar(x)&&x>0);

p.addParameter('sosBelowHz',    2, @(x)isscalar(x)&&x>0);
p.addParameter('bwFrac',     0.30, @(x)isscalar(x)&&x>0);
p.addParameter('bwMin',      0.06, @(x)isscalar(x)&&x>0);
p.addParameter('sosOrder',      4, @(x)isscalar(x)&&x>=1);
p.addParameter('sosEdgeCycles', 3, @(x)isscalar(x)&&x>=0);
p.addParameter('winCapSec',    50, @(x)isscalar(x)&&x>0);

p.parse(varargin{:});
opt = p.Results;

winShape     = lower(string(opt.winShape));
nCyclesWin   = opt.nCyclesWin;
nShuf        = opt.nShuf;
waveCycles   = opt.waveCycles;
envLP        = opt.envLowpassHz;
targetFs     = opt.targetFs;
minShiftSec  = opt.minShiftSec;
chunkSize    = opt.chunkSize;

sosBelowHz   = opt.sosBelowHz;
bwFrac       = opt.bwFrac;
bwMin        = opt.bwMin;
sosOrder     = opt.sosOrder;
sosEdgeCycles= opt.sosEdgeCycles;
winCapSec    = opt.winCapSec;

% -------------------- inputs / mode selection --------------------
useRespMode = isempty(macRaw);
disp('inputs parsed and establishing mode')
if ~isfield(rawDat,'data') || isempty(rawDat.data)
    error('rawDat.data is required.');
end

xData = double(rawDat.data(:).'); % 1 x T

if useRespMode
    if ~isfield(rawDat,'rsp') || isempty(rawDat.rsp)
        error('Resp mode requires rawDat.rsp.');
    end
    xPhaseSrc = double(rawDat.rsp(:).');   % respiration phase source
    xPhaseSrc = detrend(xPhaseSrc);        % per spec
    xAmpSrc   = xData;                     % amplitude from rawDat.data
    phaseSource = "rawDat.rsp";
    ampSource   = "rawDat.data";
else
    if ~isfield(macRaw,'data') || isempty(macRaw.data)
        error('Legacy mode requires macRaw.data.');
    end
    xPhaseSrc = xData;                     % EEG phase source
    xAmpSrc   = double(macRaw.data(:).');  % OB amplitude source
    phaseSource = "rawDat.data";
    ampSource   = "macRaw.data";
end

if numel(xPhaseSrc) ~= numel(xAmpSrc)
    error('Phase source and amplitude source must have the same length (shared time base).');
end
T = numel(xPhaseSrc);

PACfrex = PACfrex(:);
nFrex   = numel(PACfrex);

[nBreaths, nTpl] = size(keyBreathIDX);

% Fill NaNs gently (filters hate NaNs)
if any(~isfinite(xPhaseSrc))
    xPhaseSrc = fillmissing(xPhaseSrc,'linear','EndValues','nearest');
end
if any(~isfinite(xAmpSrc))
    xAmpSrc = fillmissing(xAmpSrc,'linear','EndValues','nearest');
end
disp('mode established and source vectors set')
% -------------------- amplitude envelope (bandpass -> Hilbert -> lowpass -> resample) --------------------
bp = bpHz(:).';
if numel(bp) ~= 2 || bp(1) <= 0 || bp(2) <= bp(1) || bp(2) >= fs/2
    error('bpHz must be [low high] with 0<low<high<Nyquist.');
end

bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 8, ...
    'HalfPowerFrequency1', bp(1), ...
    'HalfPowerFrequency2', bp(2), ...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');

xBP = filtfilt(bpFilt, xAmpSrc);
A   = abs(hilbert(xBP));

lpFilt = designfilt('lowpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency', min(envLP, fs/2-1), ...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');
A_lp = filtfilt(lpFilt, A);

% Resample amplitude envelope to targetFs
if abs(targetFs - fs) > 1e-12
    A_ds = resample(A_lp, targetFs, fs);
else
    A_ds = A_lp;
end
A_ds = double(A_ds(:).'); % 1 x N

% -------------------- resample phase source once (for SOS filtering) --------------------
if abs(targetFs - fs) > 1e-12
    phaseSrc_ds = resample(xPhaseSrc, targetFs, fs);
else
    phaseSrc_ds = xPhaseSrc;
end
phaseSrc_ds = double(phaseSrc_ds(:).'); % 1 x N

% Enforce same length (resample should match, but keep robust)
N = min(numel(A_ds), numel(phaseSrc_ds));
A_ds        = A_ds(1:N);
phaseSrc_ds = phaseSrc_ds(1:N);

% -------------------- adjust keyBreathIDX to targetFs grid --------------------
keyIDX_ds = nan(size(keyBreathIDX));
m = isfinite(keyBreathIDX);
keyIDX_ds(m) = round((keyBreathIDX(m)-1) * (targetFs/fs)) + 1;

bad = m & (keyIDX_ds < 1 | keyIDX_ds > N);
keyIDX_ds(bad) = NaN;

% -------------------- precompute GLOBAL null shifts --------------------
minShiftSamp = round(minShiftSec * targetFs);
if N <= 2*minShiftSamp + 10
    error('Recording too short for minShiftSec=%g at targetFs=%g (N=%d).', minShiftSec, targetFs, N);
end
rng(opt.rngSeed);
shiftsAll = randi([minShiftSamp, N-minShiftSamp], nShuf, nFrex);

% -------------------- allocate output --------------------
pacOut = nan(nBreaths, nTpl, nFrex, 9, 'single');
methodUsed = strings(nFrex,1);

% ---- global (per-frequency) metrics ----
pacGlobal = nan(nFrex, 10, 'single');      % [nFrex x 10] global versions of your 10 metrics
mzGlobal  = complex(nan(nFrex,1));         % raw complex sum vector (unnormalized)
mz0Global = complex(nan(nFrex,1));         % debiased complex sum vector (unnormalized)

% optional: how many windows contributed per frequency
nWinGlobal = zeros(nFrex,1);
disp('everything ready to begin frequency loop')
% -------------------- main loop over phase freqs --------------------
for fi = 1:nFrex
    f = PACfrex(fi);
    disp(fi)
    % decide method (hard hybrid)
    useSOS = (f < sosBelowHz);
    methodUsed(fi) = ternary(useSOS, "sos", "wavelet");

    % ----- build unit phasor U(t) at targetFs -----
    waveEdge_ds = 0;
    filtEdge_ds = 0;

    if useSOS
        % Bandwidth rule (Hz)
        bw = max(bwMin, bwFrac * f);
        fLo = f - bw/2;
        fHi = f + bw/2;

        % clamp to valid range
        fLo = max(fLo, 0.001);                   % avoid 0
        fHi = min(fHi, targetFs/2 - 0.001);      % avoid Nyquist

        if fHi <= fLo
            continue; % leave NaNs
        end

        bpFilt = designfilt('bandpassiir', ...
                    'FilterOrder', 8, ...
                    'HalfPowerFrequency1', fLo, ...
                    'HalfPowerFrequency2', fHi, ...
                    'SampleRate', targetFs, ...
                    'DesignMethod', 'butter');

        ybp = filtfilt(bpFilt, phaseSrc_ds);

        z   = hilbert(ybp);
        U   = exp(1i * angle(z));

        % conservative edge clearance: sosEdgeCycles cycles of f
        filtEdge_ds = ceil((sosEdgeCycles / max(f,0.001)) * targetFs);

    else
        % Morlet wavelet phase at original fs, then resample complex series
        s = waveCycles / (2*pi*f);
        t = (-3*s : 1/fs : 3*s);
        wv = exp(2*1i*pi*f*t) .* exp(-(t.^2)/(2*s^2));
        wv = wv ./ sqrt(sum(abs(wv).^2));

        nWave = numel(wv);
        nConv = T + nWave - 1;
        halfW = floor(nWave/2);

        Xf = fft(xPhaseSrc, nConv);
        Wf = fft(wv,        nConv);
        convx = ifft(Wf .* Xf);
        convx = convx(halfW+1 : halfW+T);

        if abs(targetFs - fs) > 1e-12
            convx_ds = resample(convx, targetFs, fs);
        else
            convx_ds = convx;
        end
        convx_ds = convx_ds(:).';
        convx_ds = convx_ds(1:N);

        U = exp(1i * angle(convx_ds));

        waveEdge_ds = ceil((nWave/2) * (targetFs/fs));
    end

    % ----- window definition (cycles -> samples @ targetFs), cap at winCapSec -----
    winLen = max(3, round((nCyclesWin / f) * targetFs));
    winCap = max(3, round(winCapSec * targetFs));
    winLen = min(winLen, winCap);
    if mod(winLen,2)==0, winLen = winLen + 1; end
    halfWin = floor(winLen/2);

    % edge clearance
    edgeClear = max([halfWin, waveEdge_ds, filtEdge_ds]);

    % window weights
    switch winShape
        case "gaussian"
            w = gausswin(winLen, 2.5);
        case "boxcar"
            w = ones(winLen,1);
        otherwise
            error('winShape must be ''gaussian'' or ''boxcar''.');
    end
    w = w / sum(w);
    w = single(w);

    % ----- index bookkeeping -----
    idxLin = keyIDX_ds(:);
    ok     = isfinite(idxLin);
    idxV   = idxLin(ok);
    posV   = find(ok);

    okEdge = (idxV > edgeClear) & (idxV <= (N - edgeClear));
    idxE   = idxV(okEdge);
    posE   = posV(okEdge);

    % linear outputs (then reshape)
    out1 = nan(nBreaths*nTpl,1,'single'); % C_A
    out2 = nan(nBreaths*nTpl,1,'single'); % C_corr
    out3 = nan(nBreaths*nTpl,1,'single'); % zA_global
    out4 = nan(nBreaths*nTpl,1,'single'); % zC_global
    out5 = nan(nBreaths*nTpl,1,'single'); % pref
    out6 = nan(nBreaths*nTpl,1,'single'); % zA_breath
    out7 = nan(nBreaths*nTpl,1,'single'); % zC_breath
    out8 = nan(nBreaths*nTpl,1,'single'); % C_debias
    out9 = nan(nBreaths*nTpl,1,'single'); % zDebias_global
    out10 = nan(nBreaths*nTpl,1,'single'); % zDebias_breath

    if isempty(idxE)
        pacOut(:,:,fi,1) = reshape(out1, nBreaths, nTpl);
        pacOut(:,:,fi,2) = reshape(out2, nBreaths, nTpl);
        pacOut(:,:,fi,3) = reshape(out3, nBreaths, nTpl);
        pacOut(:,:,fi,4) = reshape(out4, nBreaths, nTpl);
        pacOut(:,:,fi,5) = reshape(out5, nBreaths, nTpl);
        pacOut(:,:,fi,6) = reshape(out6, nBreaths, nTpl);
        pacOut(:,:,fi,7) = reshape(out7, nBreaths, nTpl);
        pacOut(:,:,fi,8) = reshape(out8, nBreaths, nTpl);
        pacOut(:,:,fi,9) = reshape(out9, nBreaths, nTpl);
        pacOut(:,:,fi,10) = reshape(out10, nBreaths, nTpl);
        continue
    end

    % NEW: precompute BREATH-wise shifts for this frequency
    shiftsBreath = randi([minShiftSamp, N-minShiftSamp], nShuf, nBreaths);

    % ---- GLOBAL accumulators for this frequency (observed + nulls) ----
    sumMz_obs     = complex(0,0);   % sum over windows of mz = <A*u>_w
    sumMz0_obs    = complex(0,0);   % sum over windows of mz0 = <(A-<A>)*u>_w (debiased)
    sumDenA_obs   = 0;             % sum over windows of <A>_w
    sumA2_obs     = 0;             % sum over windows of <A^2>_w
    sumDenDeb2_obs= 0;             % sum over windows of <(A-<A>)^2>_w  (i.e., denomDeb^2)
    nWin_obs      = 0;             % number of windows that contributed
    
    % null accumulators: numerator only (denominators don't depend on the shuffle)
    sumMz_shG   = complex(zeros(nShuf,1));   % GLOBAL shift: sum of mzSh across windows
    sumMz0_shG  = complex(zeros(nShuf,1));   % GLOBAL shift: sum of mz0Sh across windows
    
    sumMz_shB   = complex(zeros(nShuf,1));   % BREATH-wise shift: sum of mzShB across windows
    sumMz0_shB  = complex(zeros(nShuf,1));   % BREATH-wise shift: sum of mz0ShB across windows

    % chunked compute
    nE = numel(idxE);
    nChunks = ceil(nE / chunkSize);

    A_ds_local = single(A_ds(:));   % N x 1
    U_local    = U(:);              % N x 1 complex double
    w_local    = w(:);              % winLen x 1 single
    offs = int32((-halfWin:halfWin));

    for ck = 1:nChunks
        a = (ck-1)*chunkSize + 1;
        b = min(ck*chunkSize, nE);
        idxChunk = idxE(a:b);
        posChunk = posE(a:b);
        nC = numel(idxChunk);

        idxMat = int32(idxChunk) + offs;
        idxMat = max(1, min(int32(N), idxMat));

        % map linear positions -> breath indices (for breath-wise shuffles)
        
        breathIdxChunk = int32(mod(double(posChunk)-1, double(nBreaths)) + 1);
        Awin = A_ds_local(idxMat);             % [nC x winLen] single
        mA   = Awin * w_local;                 % [nC x 1] single (weighted mean)
        mA2  = (Awin.^2) * w_local;            % [nC x 1] single
        denomA    = mA;
        denomCorr = sqrt(mA2);

        Uwin = U_local(idxMat);                % [nC x winLen] complex double

        % observed (standard)
        mz    = (double(Awin) .* Uwin) * double(w_local);
        absMz = abs(mz);
        C_A    = absMz ./ max(double(denomA), eps);
        C_corr = absMz ./ max(double(denomCorr), eps);
        pref   = angle(mz);

        % observed (DEBIASED): A0 = A - <A>_w
        A0      = double(Awin) - double(mA);                 % [nC x winLen] (implicit expansion)
        denomDeb= sqrt( (A0.^2) * double(w_local) );         % RMS of debiased amplitude in window
        mz0     = (A0 .* Uwin) * double(w_local);
        C_deb   = abs(mz0) ./ max(denomDeb, eps);

        % ---- accumulate OBSERVED global sums over windows ----
        sumMz_obs      = sumMz_obs + sum(mz);
        sumMz0_obs     = sumMz0_obs + sum(mz0);
        
        sumDenA_obs    = sumDenA_obs + sum(double(mA));          % <A>_w
        sumA2_obs      = sumA2_obs   + sum(double(mA2));         % <A^2>_w
        sumDenDeb2_obs = sumDenDeb2_obs + sum(double(denomDeb).^2); % <(A-<A>)^2>_w
        
        nWin_obs       = nWin_obs + nC;



        % null accumulators (Welford)
        % GLOBAL null for C_A / C_corr / C_deb
        muA = zeros(nC,1); m2A = zeros(nC,1);
        muC = zeros(nC,1); m2C = zeros(nC,1);
        muD = zeros(nC,1); m2D = zeros(nC,1);

        % BREATH-wise null for C_A / C_corr
        muAb = zeros(nC,1); m2Ab = zeros(nC,1);
        muCb = zeros(nC,1); m2Cb = zeros(nC,1);
        % BREATH-wise null for C_A / C_corr
        muAb = zeros(nC,1); m2Ab = zeros(nC,1);
        muCb = zeros(nC,1); m2Cb = zeros(nC,1);
        % BREATH-wise null for C_deb (NEW)
        muDb = zeros(nC,1); m2Db = zeros(nC,1);

        for ss = 1:nShuf
            % -------- GLOBAL circular shift (existing) --------
            sh = shiftsAll(ss, fi);
            idxMatSh = mod(double(idxMat)-1 + sh, N) + 1;     % [nC x winLen]
            UwinSh   = U_local(idxMatSh);

            mzSh   = (double(Awin) .* UwinSh) * double(w_local);
            absMzS = abs(mzSh);

            CA_s    = absMzS ./ max(double(denomA), eps);
            Ccorr_s = absMzS ./ max(double(denomCorr), eps);

            % debiased under GLOBAL shift
            mz0Sh   = (A0 .* UwinSh) * double(w_local);
            Cdeb_s  = abs(mz0Sh) ./ max(denomDeb, eps);

            % ---- accumulate GLOBAL-shift numerators across windows ----
            sumMz_shG(ss)  = sumMz_shG(ss)  + sum(mzSh);
            sumMz0_shG(ss) = sumMz0_shG(ss) + sum(mz0Sh);

            % Welford update (GLOBAL)
            k  = ss;
            d  = CA_s    - muA;  muA = muA + d./k;   m2A = m2A + d .* (CA_s    - muA);
            d2 = Ccorr_s - muC;  muC = muC + d2./k;  m2C = m2C + d2.* (Ccorr_s - muC);
            d3 = Cdeb_s  - muD;  muD = muD + d3./k;  m2D = m2D + d3.* (Cdeb_s  - muD);

            % -------- BREATH-wise circular shift (new) --------
            shVec = double(shiftsBreath(ss, breathIdxChunk)); % [nC x 1]
            idxMatShB = zeros(size(idxMat), 'double');        % [nC x winLen]
            for rr = 1:nC
                % idxMat(rr,:) is int32 indices; convert to double for mod arithmetic
                idxMatShB(rr,:) = mod(double(idxMat(rr,:)) - 1 + shVec(rr), N) + 1;
            end
            UwinShB   = U_local(idxMatShB);

            mzShB   = (double(Awin) .* UwinShB) * double(w_local);
            absMzSB = abs(mzShB);

            CA_sb    = absMzSB ./ max(double(denomA), eps);
            Ccorr_sb = absMzSB ./ max(double(denomCorr), eps);

            % debiased under BREATH-wise shift (NEW)
            mz0ShB  = (A0 .* UwinShB) * double(w_local);
            Cdeb_sb = abs(mz0ShB) ./ max(denomDeb, eps);

            % ---- accumulate BREATH-wise-shift numerators across windows ----
            sumMz_shB(ss)  = sumMz_shB(ss)  + sum(mzShB);
            sumMz0_shB(ss) = sumMz0_shB(ss) + sum(mz0ShB);

            % Welford update (BREATH-wise)
            d  = CA_sb    - muAb; muAb = muAb + d./k;   m2Ab = m2Ab + d .* (CA_sb    - muAb);
            d2 = Ccorr_sb - muCb; muCb = muCb + d2./k;  m2Cb = m2Cb + d2.* (Ccorr_sb - muCb);
            d3 = Cdeb_sb - muDb;  muDb = muDb + d3./k;  m2Db = m2Db + d3 .* (Cdeb_sb - muDb);
        end

        % finalize null SDs
        sdA  = sqrt(m2A  ./ max(nShuf-1,1));
        sdC  = sqrt(m2C  ./ max(nShuf-1,1));
        sdD  = sqrt(m2D  ./ max(nShuf-1,1));

        sdAb = sqrt(m2Ab ./ max(nShuf-1,1));
        sdCb = sqrt(m2Cb ./ max(nShuf-1,1));
        sdDb = sqrt(m2Db ./ max(nShuf-1,1));

        % z-scores
        zA  = (C_A    - muA)  ./ max(sdA,  eps);
        zC  = (C_corr - muC)  ./ max(sdC,  eps);
        zDb = (C_deb  - muD)  ./ max(sdD,  eps);

        zAb = (C_A    - muAb) ./ max(sdAb, eps);
        zCb = (C_corr - muCb) ./ max(sdCb, eps);
        zDbb = (C_deb - muDb) ./ max(sdDb, eps);

        % write chunk
        out1(posChunk) = single(C_A);
        out2(posChunk) = single(C_corr);
        out3(posChunk) = single(zA);
        out4(posChunk) = single(zC);
        out5(posChunk) = single(pref);

        out6(posChunk) = single(zAb);
        out7(posChunk) = single(zCb);
        out8(posChunk) = single(C_deb);
        out9(posChunk) = single(zDb);
        out10(posChunk) = single(zDbb);
    end

    % ============================
    % GLOBAL (per-frequency) metrics
    % ============================
    denA  = max(sumDenA_obs, eps);
    denC  = max(sqrt(sumA2_obs), eps);
    denDb = max(sqrt(sumDenDeb2_obs), eps);
    
    % observed global metrics
    C_A_g    = abs(sumMz_obs)  / denA;
    C_corr_g = abs(sumMz_obs)  / denC;
    pref_g   = angle(sumMz_obs);
    
    C_deb_g  = abs(sumMz0_obs) / denDb;
    
    % null distributions (global shift)
    CA_nullG    = abs(sumMz_shG)  / denA;
    Ccorr_nullG = abs(sumMz_shG)  / denC;
    Cdeb_nullG  = abs(sumMz0_shG) / denDb;
    
    mu_CA_G = mean(CA_nullG,    'omitnan');  sd_CA_G = std(CA_nullG,    0, 'omitnan');
    mu_CC_G = mean(Ccorr_nullG, 'omitnan');  sd_CC_G = std(Ccorr_nullG, 0, 'omitnan');
    mu_CD_G = mean(Cdeb_nullG,  'omitnan');  sd_CD_G = std(Cdeb_nullG,  0, 'omitnan');
    
    zA_g   = (C_A_g    - mu_CA_G) ./ max(sd_CA_G, eps);
    zC_g   = (C_corr_g - mu_CC_G) ./ max(sd_CC_G, eps);
    zDeb_g = (C_deb_g  - mu_CD_G) ./ max(sd_CD_G, eps);
    
    % null distributions (breath-wise shift)
    CA_nullB    = abs(sumMz_shB)  / denA;
    Ccorr_nullB = abs(sumMz_shB)  / denC;
    Cdeb_nullB  = abs(sumMz0_shB) / denDb;
    
    mu_CA_B = mean(CA_nullB,    'omitnan');  sd_CA_B = std(CA_nullB,    0, 'omitnan');
    mu_CC_B = mean(Ccorr_nullB, 'omitnan');  sd_CC_B = std(Ccorr_nullB, 0, 'omitnan');
    mu_CD_B = mean(Cdeb_nullB,  'omitnan');  sd_CD_B = std(Cdeb_nullB,  0, 'omitnan');
    
    zAb_g   = (C_A_g    - mu_CA_B) ./ max(sd_CA_B, eps);
    zCb_g   = (C_corr_g - mu_CC_B) ./ max(sd_CC_B, eps);
    zDebb_g = (C_deb_g  - mu_CD_B) ./ max(sd_CD_B, eps);
    
    % pack into pacGlobal in SAME order as your per-window metrics
    pacGlobal(fi,:) = single([ ...
        C_A_g, C_corr_g, zA_g, zC_g, pref_g, zAb_g, zCb_g, C_deb_g, zDeb_g, zDebb_g ]);
    
    % store complex resultants too (useful for your “directionality cancels” story)
    mzGlobal(fi)  = sumMz_obs;
    mz0Global(fi) = sumMz0_obs;
    nWinGlobal(fi)= nWin_obs;

    pacOut(:,:,fi,1) = reshape(out1, nBreaths, nTpl);
    pacOut(:,:,fi,2) = reshape(out2, nBreaths, nTpl);
    pacOut(:,:,fi,3) = reshape(out3, nBreaths, nTpl);
    pacOut(:,:,fi,4) = reshape(out4, nBreaths, nTpl);
    pacOut(:,:,fi,5) = reshape(out5, nBreaths, nTpl);
    pacOut(:,:,fi,6) = reshape(out6, nBreaths, nTpl);
    pacOut(:,:,fi,7) = reshape(out7, nBreaths, nTpl);
    pacOut(:,:,fi,8) = reshape(out8, nBreaths, nTpl);
    pacOut(:,:,fi,9) = reshape(out9, nBreaths, nTpl);
    pacOut(:,:,fi,10) = reshape(out10, nBreaths, nTpl);
end
disp('packaging for exit')
% -------------------- meta --------------------
meta = struct();
meta.mode            = ternary(useRespMode, "respMode", "legacyMode");
meta.phaseSource     = char(phaseSource);
meta.ampSource       = char(ampSource);

meta.gamMed          = gamMed;
meta.bpHz            = bpHz;
meta.PACfrex         = PACfrex;

meta.fs              = fs;
meta.targetFs        = targetFs;
meta.resampleRatio   = targetFs / fs;

meta.nCyclesWin      = nCyclesWin;
meta.winCapSec       = winCapSec;
meta.winShape        = char(winShape);

meta.nShuf           = nShuf;
meta.waveCycles      = waveCycles;
meta.envLowpassHz    = envLP;
meta.minShiftSec     = minShiftSec;

meta.sosBelowHz      = sosBelowHz;
meta.bwFrac          = bwFrac;
meta.bwMin           = bwMin;
meta.sosOrder        = sosOrder;
meta.sosEdgeCycles   = sosEdgeCycles;

meta.keyBreathIDX_ds = keyIDX_ds;
meta.phaseMethodUsed = methodUsed;

meta.global = struct();
meta.global.metrics = pacGlobal;   % [nFrex x 10]
meta.global.labels  = { ...
    'C_A', 'C_corr', 'zA_globalShift', 'zC_globalShift', 'pref', ...
    'zA_breathShift', 'zC_breathShift', 'C_debias', 'zDebias_globalShift', 'zDebias_breathShift'};

meta.global.mz   = mzGlobal;       % complex unnormalized sum
meta.global.mz0  = mz0Global;      % complex unnormalized sum (debiased)
meta.global.nWin = nWinGlobal;     % how many windows contributed

% bookkeeping for new outputs
meta.pacOut_nMetrics = 9;
meta.pacOut_labels = { ...
    'C_A', 'C_corr', 'zC_A_globalShift', 'zC_corr_globalShift', 'prefPhase', ...
    'zC_A_breathShift', 'zC_corr_breathShift', 'C_debiased', 'zC_debiased_globalShift'};
disp('end of function')
end

% =========================
% Helpers
% =========================
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end