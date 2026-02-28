function out = helper_breathISPC(rawDat, varargin)
% helper_breathISPC
% Breath-by-breath ISPC (PLV) between EEG and respiration across freqs 0.11–2 Hz,
% with TWO nulls:
%   (1) full circular-shift null (existing): circularly shift respiration phasor b(t)
%   (2) block-wise permutation null (new): split b(t) into 10 s blocks, randomly permute
%       block order, then randomly circular shift the resulting timeseries
%
% Also computes GLOBAL (session-level) ISPC per frequency, plus z and p-values
% against BOTH nulls.
%
% Requires:
%   rawDat.data  [1 x T] EEG channel timeseries (single channel)
%   rawDat.rsp   [1 x T] respiration timeseries
%   rawDat.fs    sampling rate (Hz)
%   rawDat.behDat.finalOnset  [nBreath x 1] onset sample indices (in rawDat.fs samples)
%   rawDat.behDat.length      [nBreath x 1] breath length in seconds
%
% Optional:
%   rawDat.behDat.use [nBreath x 1] include mask (1=use)
%
% Outputs (breath-wise; inflated back to nBreath0):
%   out.ispcObs         [nBreath0 x nFrex]
%   out.ispcNullM       [nBreath0 x nFrex]  (circular-shift)
%   out.ispcNullSD      [nBreath0 x nFrex]
%   out.ispcZ           [nBreath0 x nFrex]
%   out.ispcP           [nBreath0 x nFrex]
%
%   out.ispcNullM_blk   [nBreath0 x nFrex]  (block-wise permute+shift)
%   out.ispcNullSD_blk  [nBreath0 x nFrex]
%   out.ispcZ_blk       [nBreath0 x nFrex]
%   out.ispcP_blk       [nBreath0 x nFrex]
%
%   out.powObs          [nBreath0 x nFrex] mean narrowband EEG power within breath
%
% Global (session-level; per frequency):
%   out.global.ispcObs        [nFrex x 1]  (|mean phase-diff over all valid breath samples|)
%   out.global.prefPhase      [nFrex x 1]  (angle of global mean phase-diff)
%   out.global.ispcZ          [nFrex x 1]  (vs circular-shift null)
%   out.global.ispcP          [nFrex x 1]
%   out.global.ispcZ_blk      [nFrex x 1]  (vs block-wise null)
%   out.global.ispcP_blk      [nFrex x 1]
%   out.global.nullM / nullSD (for each null) also provided
%
% Other:
%   out.frex, out.fsUsed, out.breathOnOff, out.validBreath, etc.

% -------------------------
% Parse params
% -------------------------
p = inputParser;
p.addParameter('frex', logspace(log10(0.11), log10(2), 40)', @(x)isnumeric(x)&&isvector(x));
p.addParameter('fsOut', 20, @(x)isnumeric(x)&&isscalar(x)&&x>0); % downsample target
p.addParameter('bwFrac', 0.30, @(x)isnumeric(x)&&isscalar(x)&&x>0); % bandwidth = max(bwMin, bwFrac*f)
p.addParameter('bwMin', 0.06, @(x)isnumeric(x)&&isscalar(x)&&x>0);  % Hz (important for low f)
p.addParameter('lpOrder', 4, @(x)isnumeric(x)&&isscalar(x)&&x>=1);   % (kept for compat; not used below)
p.addParameter('nShuf', 200, @(x)isnumeric(x)&&isscalar(x)&&x>=10);
p.addParameter('rngSeed', 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('respPhaseMode', 'match', @(s)ischar(s)||isstring(s)); % 'match' or 'broad'
p.addParameter('respBand', [0.05 3], @(x)isnumeric(x)&&numel(x)==2&&x(1)>0&&x(2)>x(1));
p.addParameter('minBreathSamples', 5, @(x)isnumeric(x)&&isscalar(x)&&x>=1); % after resampling

% NEW: block-wise null controls
p.addParameter('blockSec', 10, @(x)isnumeric(x)&&isscalar(x)&&x>0);   % block length (sec)
p.parse(varargin{:});
P = p.Results;

frex  = P.frex(:);
nFrex = numel(frex);

% -------------------------
% Pull signals, basic checks
% -------------------------
x0 = double(rawDat.data(:)');
y0 = double(rawDat.rsp(:)');
if numel(x0) ~= numel(y0)
    error('rawDat.data and rawDat.rsp must be the same length.');
end
fs0 = double(rawDat.fs);
T0  = numel(x0);

% Fill NaNs gently (filters hate NaNs)
if any(~isfinite(x0)), x0 = fillmissing(x0,'linear','EndValues','nearest'); end
if any(~isfinite(y0)), y0 = fillmissing(y0,'linear','EndValues','nearest'); end

x0 = x0 - mean(x0);
y0 = y0 - mean(y0);

% -------------------------
% Downsample (recommended for 0.11–2 Hz work)
% -------------------------
fsUsed = min(fs0, P.fsOut);
if fsUsed < fs0
    x = resample(x0, fsUsed, fs0);
    y = resample(y0, fsUsed, fs0);
else
    x = x0;
    y = y0;
end
x = x(:)'; y = y(:)';  % row vectors
T = numel(x);

% -------------------------
% Breath windows (convert to fsUsed indices)
% -------------------------
on0    = rawDat.behDat.finalOnset(:);
lenSec = rawDat.behDat.length(:);
if numel(on0) ~= numel(lenSec)
    error('behDat.finalOnset and behDat.length must have same length.');
end
nBreath0 = numel(on0);

% Optional use-mask
if isfield(rawDat.behDat,'use') && ~isempty(rawDat.behDat.use)
    useMask = rawDat.behDat.use(:)==1;
else
    useMask = true(nBreath0,1);
end

addVal = 25;
on  = round((on0-1) * (fsUsed/fs0)) + 1 - round(fsUsed*addVal);
dur = round(lenSec * fsUsed);
off = on + dur - 1 + round(fsUsed*addVal*2);

% Clamp & keep valid breaths
valid = useMask & isfinite(on) & isfinite(off) & dur>=P.minBreathSamples & on>=1 & off<=T;
on  = on(valid);
off = off(valid);
dur = dur(valid);
nBreath = numel(on);
if nBreath==0
    error('No valid breaths after applying useMask/bounds/minBreathSamples.');
end

breathOnOff = [on off];

% -------------------------
% Resp broad phase (optional mode)
% -------------------------
respMode = lower(string(P.respPhaseMode));
phiRespBroad = [];
if respMode == "broad"
    phiRespBroad = angle(hilbert(bandpass_iir(y, fsUsed, P.respBand(1), P.respBand(2), 4)));
end

% -------------------------
% Shuffle plans
% -------------------------
rng(P.rngSeed);

% circular shift null (existing)
shiftsCirc = randi([1 T-1], P.nShuf, 1); % avoid 0 shift

% block-wise permute + circular shift null (new)
blockLen = max(1, round(P.blockSec * fsUsed));
nBlocks  = floor(T / blockLen);
if nBlocks < 2
    warning('Block-wise null disabled: T=%d samples, blockLen=%d => nBlocks=%d (<2).', T, blockLen, nBlocks);
end
permBlocks = zeros(P.nShuf, max(nBlocks,1), 'uint32');
for ss = 1:P.nShuf
    if nBlocks >= 2
        permBlocks(ss,:) = uint32(randperm(nBlocks));
    else
        permBlocks(ss,:) = uint32(1);
    end
end
shiftsBlk = randi([1 T-1], P.nShuf, 1);

% -------------------------
% Allocate outputs (valid breaths only)
% -------------------------
ispcObs     = nan(nBreath, nFrex);
ispcNullM   = nan(nBreath, nFrex);
ispcNullSD  = nan(nBreath, nFrex);
ispcZ       = nan(nBreath, nFrex);
ispcP       = nan(nBreath, nFrex);

ispcNullM_blk  = nan(nBreath, nFrex);
ispcNullSD_blk = nan(nBreath, nFrex);
ispcZ_blk      = nan(nBreath, nFrex);
ispcP_blk      = nan(nBreath, nFrex);

powObs      = nan(nBreath, nFrex);

% Global (session-level) per frequency
glob_ispcObs   = nan(nFrex,1);
glob_prefPhase = nan(nFrex,1);

glob_nullM     = nan(nFrex,1);
glob_nullSD    = nan(nFrex,1);
glob_ispcZ     = nan(nFrex,1);
glob_ispcP     = nan(nFrex,1);

glob_nullM_blk  = nan(nFrex,1);
glob_nullSD_blk = nan(nFrex,1);
glob_ispcZ_blk  = nan(nFrex,1);
glob_ispcP_blk  = nan(nFrex,1);

% -------------------------
% Main loop over frequencies
% -------------------------
for fi = 1:nFrex
    f = frex(fi);

    % Bandwidth choice (Hz)
    bw = max(P.bwMin, P.bwFrac * f);
    fLo = max(f - bw/2, 0.001);
    fHi = min(f + bw/2, fsUsed/2 - 0.001);
    if fHi <= fLo
        continue
    end

    % Bandpass design for this frequency
    bpFilt = designfilt('bandpassiir', ...
        'FilterOrder', 8, ...
        'HalfPowerFrequency1', fLo, ...
        'HalfPowerFrequency2', fHi, ...
        'SampleRate', fsUsed, ...
        'DesignMethod', 'butter');

    % EEG analytic
    xbp = filtfilt(bpFilt, x);
    zX  = hilbert(xbp);
    phiX = angle(zX);
    powX = abs(zX).^2;

    % Resp phase (matched vs broad)
    if respMode == "match"
        ybp = filtfilt(bpFilt, y);
        zY  = hilbert(ybp);
        phiY = angle(zY);
    else
        phiY = phiRespBroad;
    end

    % Unit phasors and phase-diff complex vector
    a = exp( 1i*phiX);
    b = exp(-1i*phiY);
    cObs = a .* b;

    % Breath-wise ISPC via cumulative sum
    cObsC = cObs(:);           % T x 1 complex
    cs = cumsum([0; cObsC]);   % (T+1) x 1
    sumsObs = cs(off+1) - cs(on);        % nBreath x 1 complex sum within each breath window
    nSamp   = off - on + 1;              % nBreath x 1
    mBreath = sumsObs ./ nSamp;          % nBreath x 1 complex mean phase-diff

    ispcObs(:,fi) = abs(mBreath);

    % Breath-wise mean power via cumulative sum
    csP = cumsum([0; powX(:)]);
    sumsP = csP(off+1) - csP(on);
    powObs(:,fi) = sumsP ./ nSamp;

    % GLOBAL observed (across all valid breath windows)
    sumGlob = sum(sumsObs, 'omitnan');
    nGlob   = sum(nSamp,   'omitnan');
    mGlob   = sumGlob ./ max(nGlob, eps);
    glob_ispcObs(fi)   = abs(mGlob);
    glob_prefPhase(fi) = angle(mGlob);

    % -------------------------
    % Null shuffles
    % -------------------------
    shMat_circ = nan(nBreath, P.nShuf);
    shMat_blk  = nan(nBreath, P.nShuf);
    globSh_circ = nan(P.nShuf,1);
    globSh_blk  = nan(P.nShuf,1);

    for ss = 1:P.nShuf
        % ===== (1) CIRCULAR SHIFT NULL (existing): shift b only =====
        cSh = a .* circshift(b, [0 shiftsCirc(ss)]);
        csSh = cumsum([0; cSh(:)]);
        sumsSh = csSh(off+1) - csSh(on);
        mSh = sumsSh ./ nSamp;

        shMat_circ(:,ss) = abs(mSh);
        globSh_circ(ss)  = abs(sum(sumsSh) ./ max(nGlob, eps));

        % ===== (2) BLOCK-WISE PERMUTE + SHIFT NULL (new): permute b blocks, then circshift =====
        if nBlocks >= 2
            bPerm = blockperm_then_shift(b, blockLen, permBlocks(ss,:), shiftsBlk(ss));
            cShB  = a .* bPerm;
            csB   = cumsum([0; cShB(:)]);
            sumsB = csB(off+1) - csB(on);
            mB    = sumsB ./ nSamp;

            shMat_blk(:,ss) = abs(mB);
            globSh_blk(ss)  = abs(sum(sumsB) ./ max(nGlob, eps));
        else
            shMat_blk(:,ss) = NaN;
            globSh_blk(ss)  = NaN;
        end
    end

    % ---- Breath-wise z/p vs circular null ----
    mu = mean(shMat_circ, 2, 'omitnan');
    sd = std(shMat_circ, 0, 2, 'omitnan');
    sd(sd==0 | ~isfinite(sd)) = NaN;

    ispcNullM(:,fi)  = mu;
    ispcNullSD(:,fi) = sd;
    ispcZ(:,fi)      = (ispcObs(:,fi) - mu) ./ sd;
    ispcP(:,fi)      = (sum(shMat_circ >= ispcObs(:,fi), 2, 'omitnan') + 1) ./ (P.nShuf + 1);

    % ---- Breath-wise z/p vs block-wise null ----
    muB = mean(shMat_blk, 2, 'omitnan');
    sdB = std(shMat_blk, 0, 2, 'omitnan');
    sdB(sdB==0 | ~isfinite(sdB)) = NaN;

    ispcNullM_blk(:,fi)  = muB;
    ispcNullSD_blk(:,fi) = sdB;
    ispcZ_blk(:,fi)      = (ispcObs(:,fi) - muB) ./ sdB;
    ispcP_blk(:,fi)      = (sum(shMat_blk >= ispcObs(:,fi), 2, 'omitnan') + 1) ./ (P.nShuf + 1);

    % ---- Global z/p vs circular null ----
    muG = mean(globSh_circ, 'omitnan');
    sdG = std(globSh_circ, 0, 'omitnan');
    if ~isfinite(sdG) || sdG==0, sdG = NaN; end

    glob_nullM(fi)  = muG;
    glob_nullSD(fi) = sdG;
    glob_ispcZ(fi)  = (glob_ispcObs(fi) - muG) ./ max(sdG, eps);
    glob_ispcP(fi)  = (sum(globSh_circ >= glob_ispcObs(fi), 'omitnan') + 1) ./ (P.nShuf + 1);

    % ---- Global z/p vs block-wise null ----
    muGB = mean(globSh_blk, 'omitnan');
    sdGB = std(globSh_blk, 0, 'omitnan');
    if ~isfinite(sdGB) || sdGB==0, sdGB = NaN; end

    glob_nullM_blk(fi)  = muGB;
    glob_nullSD_blk(fi) = sdGB;
    glob_ispcZ_blk(fi)  = (glob_ispcObs(fi) - muGB) ./ max(sdGB, eps);
    glob_ispcP_blk(fi)  = (sum(globSh_blk >= glob_ispcObs(fi), 'omitnan') + 1) ./ (P.nShuf + 1);
end

% -------------------------
% Pack outputs (inflate back to full nBreath0 with NaNs for invalid breaths)
% -------------------------
out = struct();

out.frex   = frex;
out.fsUsed = fsUsed;

out.validBreath  = valid(:);
out.nBreath0     = nBreath0;
out.nBreathValid = nBreath;

% Breath metadata
breathOnOff_full = nan(nBreath0, 2);
breathDurS_full  = nan(nBreath0, 1);

breathOnOff_full(valid,:) = breathOnOff;
breathDurS_full(valid)    = dur(:) / fsUsed;

out.breathOnOff = breathOnOff_full;
out.breathDurS  = breathDurS_full;

% Per-breath outputs (inflate)
ispcObs_full        = nan(nBreath0, nFrex);
ispcNullM_full      = nan(nBreath0, nFrex);
ispcNullSD_full     = nan(nBreath0, nFrex);
ispcZ_full          = nan(nBreath0, nFrex);
ispcP_full          = nan(nBreath0, nFrex);

ispcNullM_blk_full  = nan(nBreath0, nFrex);
ispcNullSD_blk_full = nan(nBreath0, nFrex);
ispcZ_blk_full      = nan(nBreath0, nFrex);
ispcP_blk_full      = nan(nBreath0, nFrex);

powObs_full         = nan(nBreath0, nFrex);

ispcObs_full(valid,:)        = ispcObs;
ispcNullM_full(valid,:)      = ispcNullM;
ispcNullSD_full(valid,:)     = ispcNullSD;
ispcZ_full(valid,:)          = ispcZ;
ispcP_full(valid,:)          = ispcP;

ispcNullM_blk_full(valid,:)  = ispcNullM_blk;
ispcNullSD_blk_full(valid,:) = ispcNullSD_blk;
ispcZ_blk_full(valid,:)      = ispcZ_blk;
ispcP_blk_full(valid,:)      = ispcP_blk;

powObs_full(valid,:)         = powObs;

out.ispcObs        = ispcObs_full;
out.ispcNullM      = ispcNullM_full;
out.ispcNullSD     = ispcNullSD_full;
out.ispcZ          = ispcZ_full;
out.ispcP          = ispcP_full;

out.ispcNullM_blk  = ispcNullM_blk_full;
out.ispcNullSD_blk = ispcNullSD_blk_full;
out.ispcZ_blk      = ispcZ_blk_full;
out.ispcP_blk      = ispcP_blk_full;

out.powObs         = powObs_full;

% Global outputs
out.global = struct();
out.global.ispcObs    = glob_ispcObs;
out.global.prefPhase  = glob_prefPhase;

out.global.nullM      = glob_nullM;
out.global.nullSD     = glob_nullSD;
out.global.ispcZ      = glob_ispcZ;
out.global.ispcP      = glob_ispcP;

out.global.nullM_blk  = glob_nullM_blk;
out.global.nullSD_blk = glob_nullSD_blk;
out.global.ispcZ_blk  = glob_ispcZ_blk;
out.global.ispcP_blk  = glob_ispcP_blk;

% Record null settings
out.null = struct();
out.null.circularShifts = shiftsCirc;
out.null.blockSec   = P.blockSec;
out.null.blockLenSamp = blockLen;
out.null.nBlocks    = nBlocks;
out.null.blockPerms = permBlocks;
out.null.blockShifts= shiftsBlk;

% -------------------------
% (your existing diagnostics block can remain below unchanged)
% -------------------------
% ============================================================
% Breath-wise diagnostics: fb, delta-f, local dominance, rank
% Uses:
%   fb = 1 ./ rawDat.behDat.length
% Computes for BOTH power and ISPC (raw + z):
%   - fMax, deltaF = fMax - fb
%   - local dominance: value(fb) - median(value in fb±0.2Hz, excluding fb bin)
%   - rank/percentile of fb bin within (fb±0.2Hz) neighborhood
% Stores in out.diag.*
% ============================================================

% --- breath-by-breath respiratory frequency (Hz), full length ---
fb = 1 ./ double(rawDat.behDat.length(:));

% respect validBreath mask
fb(~out.validBreath) = NaN;

out.diag = struct();
out.diag.fb = fb;  % [nBreath0 x 1]

% --- nearest frequency bin index for each breath ---
frexRow = out.frex(:)'; % 1 x nFrex
bFidx = nan(out.nBreath0, 1);

mFB = isfinite(fb);
if any(mFB)
    [~, bFidx(mFB)] = min(abs(frexRow - fb(mFB)), [], 2);
end

out.diag.fbIdx = bFidx; % [nBreath0 x 1]

% --- pull values at the fb bin for convenience ---
bFidxSafe = bFidx;
bFidxSafe(~isfinite(bFidxSafe)) = 1;

out.diag.fb_ispcRaw = arrayfun(@(bb,ii) ispcObs_full(bb,ii), (1:numel(bFidxSafe))', bFidxSafe(:));
out.diag.fb_ispcZ   = arrayfun(@(bb,ii) ispcZ_full(bb,ii),   (1:numel(bFidxSafe))', bFidxSafe(:));
out.diag.fb_powRaw  = arrayfun(@(bb,ii) powObs_full(bb,ii),  (1:numel(bFidxSafe))', bFidxSafe(:));

% neighborhood half-width for local metrics
winHz = 0.2;

nB = out.nBreath0;
nF = numel(out.frex); %#ok<NASGU>

% allocate outputs (full-length)
out.diag.ispc = struct();
out.diag.ispc.raw = struct();
out.diag.ispc.z   = struct();

out.diag.pow = struct();
out.diag.pow.raw = struct();

% --------------------------
% ISPC RAW: fMax, deltaF
% --------------------------
M = out.ispcObs;  % [nBreath0 x nFrex]
[ispcMaxRaw, ispcMaxIdxRaw] = max(M, [], 2, 'omitnan');

fMax_ispcRaw = nan(nB,1);
m = isfinite(ispcMaxIdxRaw);
fMax_ispcRaw(m) = out.frex(ispcMaxIdxRaw(m));

deltaF_ispcRaw = fMax_ispcRaw - fb;

out.diag.ispc.raw.fMax   = fMax_ispcRaw;
out.diag.ispc.raw.deltaF = deltaF_ispcRaw;
out.diag.ispc.raw.maxVal = ispcMaxRaw;
out.diag.ispc.raw.maxIdx = ispcMaxIdxRaw;

% --------------------------
% ISPC Z: fMax, deltaF
% --------------------------
Mz = out.ispcZ;
[ispcMaxZ, ispcMaxIdxZ] = max(Mz, [], 2, 'omitnan');

fMax_ispcZ = nan(nB,1);
m = isfinite(ispcMaxIdxZ);
fMax_ispcZ(m) = out.frex(ispcMaxIdxZ(m));

deltaF_ispcZ = fMax_ispcZ - fb;

out.diag.ispc.z.fMax   = fMax_ispcZ;
out.diag.ispc.z.deltaF = deltaF_ispcZ;
out.diag.ispc.z.maxVal = ispcMaxZ;
out.diag.ispc.z.maxIdx = ispcMaxIdxZ;

% --------------------------
% POWER RAW: fMax, deltaF
% --------------------------
Mp = out.powObs;
[powMax, powMaxIdx] = max(Mp, [], 2, 'omitnan');

fMax_pow = nan(nB,1);
m = isfinite(powMaxIdx);
fMax_pow(m) = out.frex(powMaxIdx(m));

deltaF_pow = fMax_pow - fb;

out.diag.pow.raw.fMax   = fMax_pow;
out.diag.pow.raw.deltaF = deltaF_pow;
out.diag.pow.raw.maxVal = powMax;
out.diag.pow.raw.maxIdx = powMaxIdx;
% ============================================================
% Local dominance + rank within neighborhood fb±0.2 Hz
%   dominance = value_at_fb - median(neighborhood excluding fb bin)
%   rankPct   = percentile rank of value_at_fb among neighborhood
% ============================================================

% ---- local helper: percentile rank (ties handled by averaging) ----
% rankPct in [0,1], NaN if neighborhood too small
pct_rank = @(vals, x) ( ...
    (sum(vals < x) + 0.5*sum(vals == x)) ./ numel(vals) );

% --- function to compute dominance + rank for a given matrix ---
function [valFb, dom, rankPctOut, medNbr, nNbr] = local_metrics(Min, fbIn, fbIdxIn, frexVec, winHz)
    nB_ = size(Min, 1);

    valFb      = nan(nB_, 1);
    dom        = nan(nB_, 1);
    rankPctOut = nan(nB_, 1);
    medNbr     = nan(nB_, 1);
    nNbr       = nan(nB_, 1);

    for bb = 1:nB_
        fbb = fbIn(bb);
        ii  = fbIdxIn(bb);

        if ~isfinite(fbb) || ~isfinite(ii)
            continue
        end

        % neighborhood bins within fb±winHz
        nbrMask = abs(frexVec - fbb) <= winHz;
        nbrIdx  = find(nbrMask);

        if numel(nbrIdx) < 2
            continue  % need at least fb bin + 1 neighbor
        end

        % value at fb bin
        x = Min(bb, ii);
        if ~isfinite(x)
            continue
        end
        valFb(bb) = x;

        % neighbors excluding fb bin
        nbrVals = Min(bb, nbrIdx);
        nbrVals(nbrIdx == ii) = NaN;          % exclude the fb bin itself
        nbrVals = nbrVals(isfinite(nbrVals));
        if isempty(nbrVals)
            continue
        end

        medNbr(bb) = median(nbrVals);
        dom(bb)    = x - medNbr(bb);

        % rank among neighborhood INCLUDING fb bin (finite only)
        allVals = Min(bb, nbrIdx);
        allVals = allVals(isfinite(allVals));
        if numel(allVals) < 2
            continue
        end

        rankPctOut(bb) = pct_rank(allVals, x);
        nNbr(bb)       = numel(allVals);
    end
end

frexVec = out.frex(:);

% --- ISPC raw neighborhood metrics ---
[out.diag.ispc.raw.valFb, ...
 out.diag.ispc.raw.localDom, ...
 out.diag.ispc.raw.fbRankPct, ...
 out.diag.ispc.raw.medNbr, ...
 out.diag.ispc.raw.nNbr] = local_metrics(out.ispcObs, fb, bFidx, frexVec, winHz);

% --- ISPC z neighborhood metrics ---
[out.diag.ispc.z.valFb, ...
 out.diag.ispc.z.localDom, ...
 out.diag.ispc.z.fbRankPct, ...
 out.diag.ispc.z.medNbr, ...
 out.diag.ispc.z.nNbr] = local_metrics(out.ispcZ, fb, bFidx, frexVec, winHz);

% --- Power raw neighborhood metrics ---
[out.diag.pow.raw.valFb, ...
 out.diag.pow.raw.localDom, ...
 out.diag.pow.raw.fbRankPct, ...
 out.diag.pow.raw.medNbr, ...
 out.diag.pow.raw.nNbr] = local_metrics(out.powObs, fb, bFidx, frexVec, winHz);

% ------------------------------------------------------------
% Notes:
% - deltaF is signed (positive means max is higher than fb).
% - fbRankPct is percentile of fb-bin within fb±0.2Hz neighborhood (0..1).
% - localDom compares fb-bin to median of nearby bins (excluding fb-bin).
% ------------------------------------------------------------

% ==========================================================
% Helper subfunctions
% ==========================================================
    function ybp = bandpass_iir(sig, fs, fLo, fHi, ord)
        wn = [fLo fHi] / (fs/2);
        wn(wn<=0) = eps;
        wn(wn>=1) = 0.999;
        [bBP,aBP] = butter(ord, wn, 'bandpass');
        ybp = filtfilt(bBP, aBP, sig);
    end

    function b2 = blockperm_then_shift(bIn, blkLen, perm, sh)
        % bIn: 1xT row. Permute full blocks of length blkLen, keep tail as-is,
        % then circularly shift entire result by sh samples.
        TT = numel(bIn);
        nB = floor(TT / blkLen);
        L  = nB * blkLen;

        if nB < 2
            b2 = circshift(bIn, [0 sh]);
            return
        end

        main = bIn(1:L);
        main = reshape(main, blkLen, nB);   % [blkLen x nB]
        main = main(:, double(perm));       % permute blocks
        b2   = [main(:).' bIn(L+1:end)];    % append tail unchanged
        b2   = circshift(b2, [0 sh]);
    end

end