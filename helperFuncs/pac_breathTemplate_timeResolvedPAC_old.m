function [pacOut, meta] = pac_breathTemplate_timeResolvedPAC_old(rawDat, macRaw, keyBreathIDX, gamMed, fs, bpHz, PACfrex, varargin)
% pac_breathTemplate_timeResolvedPAC
%
% Time-resolved PAC sampled at breath-template indices.
%   - Phase source: rawDat.data (scalp EEG), Morlet wavelet (6 cycles), phase = angle()
%   - Amplitude source: macRaw.data (OB), Butterworth IIR bandpass (bpHz), envelope = abs(hilbert(.))
%   - Envelope lowpass: 20 Hz (zero-phase), then downsample both phase+amp to targetFs (default 100 Hz)
%   - Window length: nCyclesWin cycles (default 10) with winShape (default 'gaussian')
%   - Null: circular shifts of the PHASOR u(t)=exp(1i*phi(t)) relative to amplitude A(t), with shifts >=10 s
%   - Outputs sampled at keyBreathIDX (breaths x 50) after downsampling index adjustment
%
% Output:
%   pacOut: [nBreaths x 50 x nFrex x 5]
%       (:,:,:,1) = coupling strength: |<A*u>| / <A>
%       (:,:,:,2) = corr-style:        |<A*u>| / sqrt(<A^2>)
%       (:,:,:,3) = z-score of (1) vs circular-shift null (mean/sd over shuffles at that point)
%       (:,:,:,4) = z-score of (2)
%       (:,:,:,5) = preferred phase:   angle(<A*u>)
%
% meta: struct with useful bookkeeping.
%
% Name-value options:
%   'nCyclesWin'   (default 10)
%   'nShuf'        (default 200)
%   'winShape'     (default 'gaussian')  % 'gaussian' or 'boxcar'
%   'waveCycles'   (default 6)
%   'envLowpassHz' (default 20)
%   'targetFs'     (default 100)
%   'minShiftSec'  (default 10)
%   'rngSeed'      (default 0)
%   'chunkSize'    (default 2000)

% -------------------- parse options --------------------
p = inputParser;
p.addParameter('nCyclesWin',   20, @(x)isscalar(x)&&x>0);
p.addParameter('nShuf',       200, @(x)isscalar(x)&&x>=20);
p.addParameter('winShape', 'gaussian', @(s)ischar(s)||isstring(s));
p.addParameter('waveCycles',    6, @(x)isscalar(x)&&x>0);
p.addParameter('envLowpassHz', 20, @(x)isscalar(x)&&x>0);
p.addParameter('targetFs',    100, @(x)isscalar(x)&&x>0);
p.addParameter('minShiftSec',  10, @(x)isscalar(x)&&x>0);
p.addParameter('rngSeed',       0, @(x)isscalar(x));
p.addParameter('chunkSize',  2000, @(x)isscalar(x)&&x>0);
p.parse(varargin{:});
opt = p.Results;

winShape   = lower(string(opt.winShape));
nCyclesWin = opt.nCyclesWin;
nShuf      = opt.nShuf;
waveCycles = opt.waveCycles;
envLP      = opt.envLowpassHz;
targetFs   = opt.targetFs;
minShiftSec= opt.minShiftSec;
chunkSize  = opt.chunkSize;

% -------------------- inputs --------------------
xEEG = double(rawDat.data(:).');   % 1 x T
xOB  = double(macRaw.data(:).');   % 1 x T
if numel(xEEG) ~= numel(xOB)
    error('rawDat.data and macRaw.data must have same length (shared time base).');
end
T = numel(xEEG);

PACfrex = PACfrex(:);
nFrex   = numel(PACfrex);

[nBreaths, nTpl] = size(keyBreathIDX);
if nTpl ~= 50
    % not required, but matches your pipeline assumption
end

% -------------------- downsample factor --------------------
ds = fs / targetFs;
if abs(ds - round(ds)) > 1e-9
    error('targetFs must evenly divide fs for simple decimation. fs=%g, targetFs=%g', fs, targetFs);
end
ds = round(ds);
fs_ds = fs / ds;

% -------------------- gamma envelope (OB) --------------------
% Bandpass (Butterworth IIR) -> Hilbert amplitude -> lowpass 20 Hz -> downsample
bp = bpHz(:).';
if numel(bp) ~= 2 || bp(1) <= 0 || bp(2) <= bp(1) || bp(2) >= fs/2
    error('bpHz must be [low high] with 0<low<high<Nyquist.');
end

bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 8, ...                 % 4th order butter applied twice via filtfilt; here explicit order
    'HalfPowerFrequency1', bp(1), ...
    'HalfPowerFrequency2', bp(2), ...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');

xBP = filtfilt(bpFilt, xOB);
A   = abs(hilbert(xBP));                  % raw amplitude envelope

lpFilt = designfilt('lowpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency', min(envLP, fs/2-1), ...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');
A_lp = filtfilt(lpFilt, A);

A_ds = A_lp(1:ds:end);                    % 1 x Tds (simple decimation after lowpass)
N    = numel(A_ds);

% -------------------- adjust keyBreathIDX to downsampled grid --------------------
% Map original sample index -> nearest downsampled sample
keyIDX_ds = nan(size(keyBreathIDX));
m = isfinite(keyBreathIDX);
keyIDX_ds(m) = round((keyBreathIDX(m)-1)/ds) + 1;

% drop anything that lands outside 1..N
bad = m & (keyIDX_ds < 1 | keyIDX_ds > N);
keyIDX_ds(bad) = NaN;

% -------------------- precompute null shifts (parfor-friendly) --------------------
minShiftSamp = round(minShiftSec * fs_ds);
if N <= 2*minShiftSamp + 10
    error('Recording too short for minShiftSec=%g at fs_ds=%g (N=%d).', minShiftSec, fs_ds, N);
end

rng(opt.rngSeed);
shiftsAll = randi([minShiftSamp, N-minShiftSamp], nShuf, nFrex);

% -------------------- allocate output --------------------
pacOut = nan(nBreaths, nTpl, nFrex, 5, 'single');

% -------------------- main loop over phase freqs --------------------
for fi = 1:nFrex
    f = PACfrex(fi);
f
    % ----- Morlet wavelet phase extraction at original fs -----
    s = waveCycles / (2*pi*f);                 % Gaussian SD (sec)
    t = (-3*s : 1/fs : 3*s);                   % sec
    wv = exp(2*1i*pi*f*t) .* exp(-(t.^2)/(2*s^2));
    wv = wv ./ sqrt(sum(abs(wv).^2));          % energy normalize

    nWave = numel(wv);
    nConv = T + nWave - 1;
    halfW = floor(nWave/2);

    Xf = fft(xEEG, nConv);
    Wf = fft(wv,   nConv);
    convx = ifft(Wf .* Xf);
    convx = convx(halfW+1 : halfW+T);          % 1 x T (analytic-ish bandpassed)

    % downsample complex time series (band-limited around f; safe for 100 Hz Nyquist=50)
    convx_ds = convx(1:ds:end);                % 1 x N
    phi = angle(convx_ds);
    U   = exp(1i*phi);                         % 1 x N unit phasor

    % ----- window definition (cycles -> samples @ fs_ds) -----
    winLen = max(3, round((nCyclesWin / f) * fs_ds)); % samples
    if mod(winLen,2)==0, winLen = winLen + 1; end
    halfWin = floor(winLen/2);

    % add extra clearance for wavelet edge effects after downsampling
    waveEdge_ds = ceil((nWave/2) / ds);
    edgeClear = max(halfWin, waveEdge_ds);

    % window weights
    switch winShape
        case "gaussian"
            w = gausswin(winLen, 2.5);
        case "boxcar"
            w = ones(winLen,1);
        otherwise
            error('winShape must be ''gaussian'' or ''boxcar''.');
    end
    w = w / sum(w);                            % normalize to sum=1
    w = single(w);

    % ----- prepare index bookkeeping for this frequency -----
    idxLin = keyIDX_ds(:);                     % [nBreaths*50 x 1], may contain NaN
    ok     = isfinite(idxLin);
    idxV   = idxLin(ok);
    posV   = find(ok);

    % enforce edge clearance (NaN outputs if not enough room)
    okEdge = (idxV > edgeClear) & (idxV <= (N - edgeClear));
    idxE   = idxV(okEdge);
    posE   = posV(okEdge);

    % linear buffers for outputs (then reshape to breaths x 50)
    out1 = nan(nBreaths*nTpl,1,'single'); % strength_A
    out2 = nan(nBreaths*nTpl,1,'single'); % strength_corr
    out3 = nan(nBreaths*nTpl,1,'single'); % z_strength_A
    out4 = nan(nBreaths*nTpl,1,'single'); % z_strength_corr
    out5 = nan(nBreaths*nTpl,1,'single'); % pref phase

    if isempty(idxE)
        pacOut(:,:,fi,1) = reshape(out1, nBreaths, nTpl);
        pacOut(:,:,fi,2) = reshape(out2, nBreaths, nTpl);
        pacOut(:,:,fi,3) = reshape(out3, nBreaths, nTpl);
        pacOut(:,:,fi,4) = reshape(out4, nBreaths, nTpl);
        pacOut(:,:,fi,5) = reshape(out5, nBreaths, nTpl);
        continue
    end

    % process in chunks to limit memory
    nE = numel(idxE);
    nChunks = ceil(nE / chunkSize);

    A_ds_local = single(A_ds(:));              % N x 1
    U_local    = U(:);                         % N x 1 complex double (keep as double for stability)
    w_local    = w(:);                         % winLen x 1 single

    for ck = 1:nChunks
        a = (ck-1)*chunkSize + 1;
        b = min(ck*chunkSize, nE);
        idxChunk = idxE(a:b);
        posChunk = posE(a:b);
        nC = numel(idxChunk);

        % window index matrix: [nC x winLen]
        offs = int32((-halfWin:halfWin));
        idxMat = int32(idxChunk) + offs;       % implicit expansion
        idxMat = max(1, min(int32(N), idxMat));% (shouldn't hit due to edgeClear, but safe)

        % pull amplitude window once
        Awin = A_ds_local(idxMat);             % [nC x winLen] single
        mA   = Awin * w_local;                 % [nC x 1]
        mA2  = (Awin.^2) * w_local;            % [nC x 1]
        denomA    = mA;
        denomCorr = sqrt(mA2);

        % observed
        Uwin = U_local(idxMat);                % [nC x winLen] complex double
        mz   = (double(Awin) .* Uwin) * double(w_local);  % [nC x 1] complex double
        absMz = abs(mz);

        C_A    = absMz ./ max(double(denomA), eps);
        C_corr = absMz ./ max(double(denomCorr), eps);
        pref   = angle(mz);

        % null: Welford online mean/sd at these nC points
        muA   = zeros(nC,1);
        m2A   = zeros(nC,1);
        muC   = zeros(nC,1);
        m2C   = zeros(nC,1);

        for ss = 1:nShuf
            sh = shiftsAll(ss, fi);

            % shift phasor relative to amplitude: U(t+sh)
            idxMatSh = mod(double(idxMat)-1 + sh, N) + 1; % [nC x winLen] double
            UwinSh   = U_local(idxMatSh);                 % complex double

            mzSh   = (double(Awin) .* UwinSh) * double(w_local);
            absMzS = abs(mzSh);

            CA_s    = absMzS ./ max(double(denomA), eps);
            Ccorr_s = absMzS ./ max(double(denomCorr), eps);

            % Welford update (vectorized)
            k = ss;
            d  = CA_s - muA;    muA = muA + d./k;    m2A = m2A + d.*(CA_s - muA);
            d2 = Ccorr_s - muC; muC = muC + d2./k;   m2C = m2C + d2.*(Ccorr_s - muC);
        end

        sdA = sqrt(m2A ./ max(nShuf-1,1));
        sdC = sqrt(m2C ./ max(nShuf-1,1));

        zA = (C_A    - muA) ./ max(sdA, eps);
        zC = (C_corr - muC) ./ max(sdC, eps);

        % write chunk into linear outputs
        out1(posChunk) = single(C_A);
        out2(posChunk) = single(C_corr);
        out3(posChunk) = single(zA);
        out4(posChunk) = single(zC);
        out5(posChunk) = single(pref);
    end

    % reshape into breaths x 50 and store
    pacOut(:,:,fi,1) = reshape(out1, nBreaths, nTpl);
    pacOut(:,:,fi,2) = reshape(out2, nBreaths, nTpl);
    pacOut(:,:,fi,3) = reshape(out3, nBreaths, nTpl);
    pacOut(:,:,fi,4) = reshape(out4, nBreaths, nTpl);
    pacOut(:,:,fi,5) = reshape(out5, nBreaths, nTpl);
end

% -------------------- meta --------------------
meta = struct();
meta.gamMed     = gamMed;
meta.bpHz       = bpHz;
meta.PACfrex    = PACfrex;
meta.fs         = fs;
meta.targetFs   = fs_ds;
meta.dsFactor   = ds;
meta.nCyclesWin = nCyclesWin;
meta.winShape   = char(winShape);
meta.nShuf      = nShuf;
meta.waveCycles = waveCycles;
meta.envLowpassHz = envLP;
meta.minShiftSec  = minShiftSec;
meta.keyBreathIDX_ds = keyIDX_ds;

end