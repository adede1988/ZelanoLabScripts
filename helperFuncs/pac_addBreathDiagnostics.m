function meta = pac_addBreathDiagnostics(meta, pacOut, behDat, varargin)
% pac_addBreathDiagnostics
%
% Augment PAC meta struct with breath-wise diagnostics analogous to helper_breathISPC:
%   - fb (Hz) = 1 ./ behDat.length
%   - nearest frequency bin index (fbIdx)
%   - value at fb bin (raw PAC, z PAC, preferred phase) at a chosen template index
%   - fMax and deltaF for raw PAC and z PAC (based on values at that template index)
%   - local dominance and percentile rank of fb-bin within fb ± winHz neighborhood
%
% Inputs
%   meta   : struct returned by pac_breathTemplate_timeResolvedPAC (must contain meta.PACfrex)
%   pacOut : [nBreath x nTpl x nFrex x 5] PAC output
%   behDat : struct with field .length [nBreath x 1] in seconds
%
% Name-value options
%   'tplIdxDiag' (default 2)   : which template column to use for frequency-wise diagnostics
%   'winHz'      (default 0.2) : neighborhood half-width for local dominance/rank
%   'validBreath'(default [])  : optional logical mask [nBreath x 1]; invalid -> NaN fb and outputs

% ------------------------- parse options -------------------------
p = inputParser;
p.addParameter('tplIdxDiag', 2, @(x)isscalar(x)&&x>=1);
p.addParameter('winHz', 0.2, @(x)isscalar(x)&&x>0);
p.addParameter('validBreath', [], @(x)islogical(x)||isempty(x));
p.parse(varargin{:});
opt = p.Results;

tplIdx = opt.tplIdxDiag;
winHz  = opt.winHz;

% ------------------------- checks -------------------------
if ~isfield(meta,'PACfrex') || isempty(meta.PACfrex)
    error('meta.PACfrex is required.');
end

frexVec = double(meta.PACfrex(:));         % [nFrex x 1]
nB  = size(pacOut,1);
nTpl= size(pacOut,2);
nF  = size(pacOut,3);

if tplIdx > nTpl
    error('tplIdxDiag=%d exceeds nTpl=%d in pacOut.', tplIdx, nTpl);
end

lenSec = double(behDat.length(:));
if numel(lenSec) ~= nB
    error('Length mismatch: numel(behDat.length)=%d but size(pacOut,1)=%d.', numel(lenSec), nB);
end

validBreath = opt.validBreath;
if isempty(validBreath)
    validBreath = true(nB,1);
else
    validBreath = validBreath(:);
    if numel(validBreath) ~= nB
        error('validBreath must be [nBreath x 1] to match pacOut first dimension.');
    end
end

% ------------------------- fb and nearest bin -------------------------
fb = 1 ./ lenSec;               % Hz
fb(~validBreath) = NaN;

fbIdx = nan(nB,1);
mFB = isfinite(fb);
if any(mFB)
    % nearest frequency bin index for each breath
    [~, fbIdx(mFB)] = min(abs(frexVec(:)' - fb(mFB)), [], 2);
end

% safe indices for extraction (we’ll NaN out after)
fbIdxSafe = fbIdx;
fbIdxSafe(~isfinite(fbIdxSafe)) = 1;

% ------------------------- extract value at fb bin (at tplIdx) -------------------------
% Work on 3D slices [nB x nTpl x nF]
pacRaw3 = pacOut(:,:,:,1);
pacZ3   = pacOut(:,:,:,3);
phi3    = pacOut(:,:,:,5);

linFb = sub2ind([nB, nTpl, nF], (1:nB)', repmat(tplIdx,nB,1), fbIdxSafe);

fb_pacRaw = pacRaw3(linFb);
fb_pacZ   = pacZ3(linFb);
fb_phase  = phi3(linFb);

% apply validity
fb_pacRaw(~mFB) = NaN;
fb_pacZ(~mFB)   = NaN;
fb_phase(~mFB)  = NaN;

% ------------------------- frequency-wise matrices for this template -------------------------
Mraw = squeeze(pacOut(:, tplIdx, :, 1));   % [nB x nF]
Mz   = squeeze(pacOut(:, tplIdx, :, 3));   % [nB x nF]

% ------------------------- fMax and deltaF (RAW) -------------------------
[maxValRaw, maxIdxRaw] = max(Mraw, [], 2, 'omitnan');   % [nB x 1]
fMaxRaw = nan(nB,1);
m = isfinite(maxIdxRaw);
fMaxRaw(m) = frexVec(maxIdxRaw(m));
deltaFRaw = fMaxRaw - fb;

% preferred phase at max RAW (same tplIdx, same freq bin maxIdxRaw)
maxIdxRawSafe = maxIdxRaw; maxIdxRawSafe(~isfinite(maxIdxRawSafe)) = 1;
linMaxRaw = sub2ind([nB, nTpl, nF], (1:nB)', repmat(tplIdx,nB,1), maxIdxRawSafe);
phaseAtMaxRaw = phi3(linMaxRaw);
phaseAtMaxRaw(~m) = NaN;

% ------------------------- fMax and deltaF (Z) -------------------------
[maxValZ, maxIdxZ] = max(Mz, [], 2, 'omitnan');         % [nB x 1]
fMaxZ = nan(nB,1);
m = isfinite(maxIdxZ);
fMaxZ(m) = frexVec(maxIdxZ(m));
deltaFZ = fMaxZ - fb;

maxIdxZSafe = maxIdxZ; maxIdxZSafe(~isfinite(maxIdxZSafe)) = 1;
linMaxZ = sub2ind([nB, nTpl, nF], (1:nB)', repmat(tplIdx,nB,1), maxIdxZSafe);
phaseAtMaxZ = phi3(linMaxZ);
phaseAtMaxZ(~m) = NaN;

% ------------------------- local dominance + rank within fb±winHz -------------------------
[pacRaw_valFb, pacRaw_localDom, pacRaw_fbRankPct, pacRaw_medNbr, pacRaw_nNbr] = ...
    local_metrics(Mraw, fb, fbIdx, frexVec, winHz);

[pacZ_valFb, pacZ_localDom, pacZ_fbRankPct, pacZ_medNbr, pacZ_nNbr] = ...
    local_metrics(Mz, fb, fbIdx, frexVec, winHz);

% ------------------------- pack into meta.diag -------------------------
if ~isfield(meta,'diag') || isempty(meta.diag)
    meta.diag = struct();
end

meta.diag.tplIdxDiag = tplIdx;
meta.diag.winHz      = winHz;

meta.diag.fb     = fb;       % [nB x 1]
meta.diag.fbIdx  = fbIdx;    % [nB x 1]

% quick “at fb” pulls (like your scratch)
meta.diag.fb_pacRaw = fb_pacRaw;
meta.diag.fb_pacZ   = fb_pacZ;
meta.diag.fb_phase  = fb_phase;

% structured PAC diagnostics
meta.diag.pac = struct();

meta.diag.pac.raw = struct();
meta.diag.pac.raw.fMax      = fMaxRaw;
meta.diag.pac.raw.deltaF    = deltaFRaw;
meta.diag.pac.raw.maxVal    = maxValRaw;
meta.diag.pac.raw.maxIdx    = maxIdxRaw;
meta.diag.pac.raw.phaseAtMax= phaseAtMaxRaw;

meta.diag.pac.raw.valFb     = pacRaw_valFb;
meta.diag.pac.raw.localDom  = pacRaw_localDom;
meta.diag.pac.raw.fbRankPct = pacRaw_fbRankPct;
meta.diag.pac.raw.medNbr    = pacRaw_medNbr;
meta.diag.pac.raw.nNbr      = pacRaw_nNbr;

meta.diag.pac.z = struct();
meta.diag.pac.z.fMax        = fMaxZ;
meta.diag.pac.z.deltaF      = deltaFZ;
meta.diag.pac.z.maxVal      = maxValZ;
meta.diag.pac.z.maxIdx      = maxIdxZ;
meta.diag.pac.z.phaseAtMax  = phaseAtMaxZ;

meta.diag.pac.z.valFb       = pacZ_valFb;
meta.diag.pac.z.localDom    = pacZ_localDom;
meta.diag.pac.z.fbRankPct   = pacZ_fbRankPct;
meta.diag.pac.z.medNbr      = pacZ_medNbr;
meta.diag.pac.z.nNbr        = pacZ_nNbr;

end

% =========================================================
% Local helper: dominance + rank in fb±winHz neighborhood
% =========================================================
function [valFb, dom, rankPctOut, medNbr, nNbr] = local_metrics(Min, fbIn, fbIdxIn, frexVec, winHz)
% Min: [nB x nF]
nB = size(Min,1);
valFb = nan(nB,1);
dom   = nan(nB,1);
rankPctOut = nan(nB,1);
medNbr = nan(nB,1);
nNbr = nan(nB,1);

pct_rank = @(vals, x) ( (sum(vals < x) + 0.5*sum(vals == x)) ./ numel(vals) );

for bb = 1:nB
    fbb = fbIn(bb);
    ii  = fbIdxIn(bb);
    if ~isfinite(fbb) || ~isfinite(ii), continue; end

    nbrMask = abs(frexVec - fbb) <= winHz;
    nbrIdx  = find(nbrMask);
    if numel(nbrIdx) < 2, continue; end

    x = Min(bb, ii);
    if ~isfinite(x), continue; end
    valFb(bb) = x;

    nbrVals = Min(bb, nbrIdx);
    nbrVals(nbrIdx == ii) = NaN;
    nbrVals = nbrVals(isfinite(nbrVals));
    if isempty(nbrVals), continue; end

    medNbr(bb) = median(nbrVals);
    dom(bb)    = x - medNbr(bb);

    allVals = Min(bb, nbrIdx);
    allVals = allVals(isfinite(allVals));
    if numel(allVals) < 2, continue; end

    rankPctOut(bb) = pct_rank(allVals, x);
    nNbr(bb) = numel(allVals);
end
end