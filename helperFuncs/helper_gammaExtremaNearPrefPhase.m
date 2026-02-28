function out = helper_gammaExtremaNearPrefPhase(rawDat, pacInfo, respFrex, varargin)
% helper_gammaExtremaNearPrefPhase
%
% Same as prior version, EXCEPT:
%   Within each breath, we:
%     1) find all samples with phase within +/- phaseWinRad of preferred phase
%     2) identify the LONGEST CONTIGUOUS (connected) run of those samples
%     3) search for max/min gamma envelope ONLY within that longest run
%
% (All other behavior unchanged.)

% ------------------------- parse options -------------------------
p = inputParser;
p.addParameter('phaseWinRad', pi/4, @(x)isscalar(x)&&x>0);
p.addParameter('padLocal', 1000, @(x)isscalar(x)&&x>=0);
p.addParameter('filterOrder', 8, @(x)isscalar(x)&&x>=2);
p.addParameter('gammaBPOrder', 8, @(x)isscalar(x)&&x>=2);
p.addParameter('gammaLPOrder', 4, @(x)isscalar(x)&&x>=1);
p.addParameter('fs', [], @(x)isempty(x)||(isscalar(x)&&x>0));
p.parse(varargin{:});
opt = p.Results;

phaseWin = opt.phaseWinRad;
padLocal = opt.padLocal;

% ------------------------- pull required fields -------------------------
if isempty(opt.fs)
    if isfield(rawDat,'fs') && ~isempty(rawDat.fs)
        fs = double(rawDat.fs);
    else
        error('Need fs: provide rawDat.fs or pass ''fs'', value.');
    end
else
    fs = double(opt.fs);
end

if ~isfield(rawDat,'data') || isempty(rawDat.data), error('rawDat.data required'); end
if ~isfield(rawDat,'rsp')  || isempty(rawDat.rsp),  error('rawDat.rsp required');  end

if ~isfield(pacInfo,'global') || ~isfield(pacInfo.global,'metrics')
    error('pacInfo.global.metrics required');
end

respFrex = double(respFrex(:));
nFrex = numel(respFrex);

metrics = pacInfo.global.metrics;
if size(metrics,1) ~= nFrex
    error('respFrex length (%d) must match size(pacInfo.global.metrics,1) (%d).', nFrex, size(metrics,1));
end
if size(metrics,2) < 5
    error('pacInfo.global.metrics must have at least 5 columns (pref phase at col 5).');
end
prefPhaseByF = double(metrics(:,5)); % radians

% gamma params from pacInfo
if ~isfield(pacInfo,'bpHz') || isempty(pacInfo.bpHz)
    error('pacInfo.bpHz required (gamma band [low high]).');
end
bpGamma = double(pacInfo.bpHz(:).');
if numel(bpGamma)~=2, error('pacInfo.bpHz must be [low high]'); end

if ~isfield(pacInfo,'envLowpassHz') || isempty(pacInfo.envLowpassHz)
    error('pacInfo.envLowpassHz required.');
end
envLP = double(pacInfo.envLowpassHz);

if ~isfield(pacInfo,'bwFrac') || isempty(pacInfo.bwFrac)
    error('pacInfo.bwFrac required for respiration band selection.');
end
bwFrac = double(pacInfo.bwFrac);

if ~isfield(pacInfo,'bwMin') || isempty(pacInfo.bwMin)
    error('pacInfo.bwMin required for respiration band selection.');
end
bwMin = double(pacInfo.bwMin);

% signals
x = double(rawDat.data(:).');
rsp = double(rawDat.rsp(:).');
T = numel(x);
if numel(rsp) ~= T
    error('rawDat.data and rawDat.rsp must have same length.');
end

% Fill NaNs gently
if any(~isfinite(x)),   x   = fillmissing(x,'linear','EndValues','nearest'); end
if any(~isfinite(rsp)), rsp = fillmissing(rsp,'linear','EndValues','nearest'); end

% ------------------------- breath windows -------------------------
on0 = double(rawDat.behDat.finalOnset(:));
lenSec = double(rawDat.behDat.length(:));
nB = numel(on0);
if numel(lenSec) ~= nB
    error('behDat.finalOnset and behDat.length must have same length.');
end

durSamp = round(lenSec * fs);
off0 = on0 + durSamp - 1;

on  = max(1, min(T, round(on0)));
off = max(1, min(T, round(off0)));

fb = 1 ./ lenSec;

% choose nearest respFrex per breath (tie -> lower)
fChooseIdx = nan(nB,1);
fChooseHz  = nan(nB,1);
for bb = 1:nB
    if ~isfinite(fb(bb)) || ~isfinite(on(bb)) || ~isfinite(off(bb)) || off(bb) <= on(bb)
        continue
    end
    d = abs(respFrex - fb(bb));
    m = min(d);
    ii = find(d==m, 1, 'first');
    fChooseIdx(bb) = ii;
    fChooseHz(bb)  = respFrex(ii);
end

% ------------------------- gamma envelope at fs -------------------------
bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', opt.gammaBPOrder, ...
    'HalfPowerFrequency1', bpGamma(1), ...
    'HalfPowerFrequency2', bpGamma(2), ...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');
xBP = filtfilt(bpFilt, x);
A   = abs(hilbert(xBP));

lpFilt = designfilt('lowpassiir', ...
    'FilterOrder', opt.gammaLPOrder, ...
    'HalfPowerFrequency', min(envLP, fs/2-1), ...
    'SampleRate', fs, ...
    'DesignMethod', 'butter');
A_lp = filtfilt(lpFilt, A);

% ------------------------- respiration phase per frequency (cache) -------------------------
phiByF = nan(nFrex, T, 'single');

for fi = 1:nFrex
    f = respFrex(fi);
    if ~isfinite(f) || f<=0, continue; end

    bw = max(bwMin, bwFrac * f);
    fLo = max(f - bw/2, 0.001);
    fHi = min(f + bw/2, fs/2 - 0.001);
    if fHi <= fLo, continue; end

    bpFiltR = designfilt('bandpassiir', ...
        'FilterOrder', opt.filterOrder, ...
        'HalfPowerFrequency1', fLo, ...
        'HalfPowerFrequency2', fHi, ...
        'SampleRate', fs, ...
        'DesignMethod', 'butter');

    ybp = filtfilt(bpFiltR, rsp);
    z   = hilbert(ybp);
    phiByF(fi,:) = single(angle(z));
end

% ------------------------- per-breath extrema near preferred phase (LONGEST CONTIGUOUS RUN) -------------------------
peakIDX_full  = nan(nB,1);
peakIDX_local = nan(nB,1);
minIDX_full   = nan(nB,1);
minIDX_local  = nan(nB,1);

for bb = 1:nB
    fi = fChooseIdx(bb);
    if ~isfinite(fi), continue; end

    o = on(bb);
    e = off(bb);
    if ~(isfinite(o) && isfinite(e) && e > o), continue; end

    phi = double(phiByF(fi, o:e));
    if all(~isfinite(phi)), continue; end

    pref = prefPhaseByF(fi);

    % samples within +/- phaseWin of pref (circular distance)
    dphi = angle(exp(1i*(phi - pref)));
    m = isfinite(dphi) & (abs(dphi) <= phaseWin);
    if ~any(m), continue; end

    % ---- find LONGEST CONTIGUOUS TRUE RUN in m ----
    % starts where diff goes 0->1, ends where diff goes 1->0
    dm = diff([false, m, false]);
    starts = find(dm == 1);
    ends   = find(dm == -1) - 1;
    runLen = ends - starts + 1;

    [~, iBest] = max(runLen);
    s0 = starts(iBest);
    e0 = ends(iBest);

    idxRel = s0:e0;  % indices within breath segment (continuous)

    % gamma segment
    Aseg = A_lp(o:e);
    vals = Aseg(idxRel);

    if all(~isfinite(vals)), continue; end

    % robust max/min ignoring NaN/Inf
    vMax = vals; vMax(~isfinite(vMax)) = -Inf;
    vMin = vals; vMin(~isfinite(vMin)) = +Inf;

    [~, imx] = max(vMax);
    [~, imn] = min(vMin);

    idxFullMax = o + idxRel(imx) - 1;
    idxFullMin = o + idxRel(imn) - 1;

    peakIDX_full(bb) = idxFullMax;
    minIDX_full(bb)  = idxFullMin;

    peakIDX_local(bb) = (idxFullMax - o) + padLocal;
    minIDX_local(bb)  = (idxFullMin - o) + padLocal;
end

% ------------------------- pack output -------------------------
out = struct();
out.peakIDX_full  = peakIDX_full;
out.peakIDX_local = peakIDX_local;
out.minIDX_full   = minIDX_full;
out.minIDX_local  = minIDX_local;

out.fChooseIdx = fChooseIdx;
out.fChooseHz  = fChooseHz;
out.fb         = fb;

out.meta = struct();
out.meta.respFrex = respFrex;
out.meta.phaseWinRad = phaseWin;
out.meta.padLocal = padLocal;
out.meta.fs = fs;
out.meta.gamma_bpHz = bpGamma;
out.meta.envLowpassHz = envLP;
out.meta.bwFrac = bwFrac;
out.meta.bwMin  = bwMin;
out.meta.prefPhaseByF = prefPhaseByF;

end