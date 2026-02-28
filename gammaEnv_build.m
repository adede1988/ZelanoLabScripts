function gammaEnv = gammaEnv_build(chanDat, chanFiles, filei, stem, W, pctl, nShuf)
% gammaEnv_build
% Wraps the provided working script into a function that returns gammaEnv.
% - Removes ALL plotting.
% - Keeps the same logic and heavy intermediate arrays (allPhase/allPow/allGam).
% - ONLY logic change: ITPC + null ITPC are breath-weighted (equal breath representation as much as possible),
%   and null snippets are resampled from the original allPhase breath×time×freq array (breaks peak↔breath link).
% - Adds a few helpful saved fields for later subsetting/aggregation.

if nargin < 5 || isempty(W),     W     = 100; end
if nargin < 6 || isempty(pctl),  pctl  = 98;  end
if nargin < 7 || isempty(nShuf), nShuf = 1000; end

% paths
[~, baseName] = fileparts(chanFiles(filei).name);
saveDirTf   = fullfile(stem,'CHANDAT_processed','tf_files');
gamFoo   = chanDat.fooof.gamma_peaks;
gamMax  = chanDat.fooof.gamma_peak_freq;

useVec   = double(chanDat.use(:)) == 1;

nBreaths = size(gamFoo,1);
frex = chanDat.tf.frex;
powFile = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, 1));
S = load(powFile, 'tfInfo');
gammaEnv = struct;

ampFrex = logspace(log10(.1), log10(15), 20);
allPow = nan([flip(size(S.tfInfo.pow)), 20]);
allPhase = allPow;
allGam = nan(nBreaths, size(S.tfInfo.pow,1));
gamFreq = nan(nBreaths, 1);

for ii = 1:nBreaths

    idx = find(isnan(gamFoo(ii,:)));
    gamFoo(ii,idx) = gamMax(ii,idx);
    curGam = median(gamFoo(ii,:));
    gamFreq(ii) = curGam;
    fi = find(min(abs(frex - curGam)) == abs(frex - curGam), 1);
    if ~isempty(fi)

        powFile = fullfile(saveDirTf, sprintf('%s_pow_fi%03d.mat', baseName, fi));
        S = load(powFile, 'tfInfo');
        powLow = double(S.tfInfo.pow);   % [nT x nBreaths]
        allGam(ii,:) = powLow(:,ii);
        powLow = log(powLow);

        % extract power and phase information
        powLow = powLow(:,ii) - mean(powLow(:,ii));
        [phase,pow] = multiphasevec3(ampFrex,...
                        powLow,S.tfInfo.fs,6, true);
        allPow(ii, :, :) = squeeze(pow)';
        allPhase(ii,:,:) = squeeze(phase)';

    end
end

gammaEnv.pow = nan([size(chanDat.targIDX), length(ampFrex)]);
gammaEnv.phase = nan([size(chanDat.targIDX), length(ampFrex)]);
gammaEnv.gamEnv = nan([size(chanDat.targIDX)]);

for ii = 1:nBreaths
    curIDX = round(chanDat.targIDX(ii,:) * (S.tfInfo.fs/ chanDat.fs));
    curIDX(isnan(curIDX)) = []; 
    if ~isempty(curIDX)
        curPow   = squeeze(allPow(ii,:,:));
        curPhase = squeeze(allPhase(ii,:,:));
        curEnv   = squeeze(allGam(ii,:));
       
        gammaEnv.pow(ii,1:length(curIDX),:)    = curPow(curIDX,:);
        gammaEnv.phase(ii,1:length(curIDX),:)  = curPhase(curIDX,:);
        gammaEnv.gamEnv(ii,1:length(curIDX)) = curEnv(curIDX);
     
    end
end

gammaEnv.gamFreq = gamFreq;
gammaEnv.frex    = ampFrex;
gammaEnv.fs      = S.tfInfo.fs;

% Helpful extra storage for later
gammaEnv.fs_raw  = chanDat.fs;
gammaEnv.W       = W;
gammaEnv.pctl    = pctl;
gammaEnv.nShuf   = nShuf;

% ============================================================
% Gamma-power up-crossings @ pctl percentile -> event snippets
% ============================================================
fs_tf = chanDat.fs / 5; % (kept as in working code)

% ---- 1) threshold (global across all use==1 breaths & timepoints) ----
mUse = useVec(:)==1;
x = allGam(mUse, :);
x = x(isfinite(x));                 % drop NaNs/Infs
thr = prctile(x, pctl);

% ---- 2) find up-crossings per breath ----
evBreath = [];
evT0     = [];

nT = size(allGam,2);
for ii = 1:nBreaths
    if ~mUse(ii), continue; end
    g = allGam(ii,:);
    if all(~isfinite(g)), continue; end

    above = g > thr;

    % rising edge: (t-1)<=thr and (t)>thr
    t0 = find( above(2:end) & ~above(1:end-1) ) + 1;

    for tt = 1:length(t0)

        t1 = t0(tt) + find( above(t0(tt):end-1) & ~above(t0(tt)+1:end),1) - 1;
        if isempty(t1)
            t1 = length(g);
        end
        [~, t1] = max(g(t0(tt):t1));
        t0(tt) = t1 + t0(tt) - 1;
    end
    % keep only events with full +/-W window in-bounds
    t0 = t0(t0 > W & t0 <= (nT - W));

    rawVals = chanDat.trial.data(ii, t0*5);
    t0(abs(rawVals)>100) = [] ;

    if ~isempty(t0)
        evBreath = [evBreath; repmat(ii, numel(t0), 1)];
        evT0     = [evT0;     t0(:)];
    end
end

fprintf('Found %d up-crossings above %g (%dth pct).\n', numel(evT0), thr, pctl);

%%%%% CHANGE: 
% evT0 = chanDat.gammaBurst.t0_idx; 
% evBreath = 1:length(evT0); 
% evBreath(isnan(evT0)) = []; 
% evT0(isnan(evT0)) = []; 
% evBreath = evBreath(:); 
% evT0 = evT0(:); 
%%%%%% END CHANGE

% ---- 3) extract phase snippets around each event ----
nE   = numel(evT0);
nF   = size(allPhase,3);            % should be 20 (ampFrex)
nOff = 2*W + 1;
offs = -W:W;

if nE == 0
    warning('No events found. Consider lowering percentile or checking allGam values.');
    phSnip = nan(0, nOff, nF);
    powSnip= nan(0, nOff, nF);
    rawSnip= nan(0, nOff);
    envSnip= nan(0, nOff);
else
    phSnip = nan(nE, nOff, nF);
    powSnip= nan(nE, nOff, nF);
    rawSnip= nan(nE, nOff);
    envSnip= nan(nE, nOff);

    for e = 1:nE
        ii  = evBreath(e);
        t0  = evT0(e);
        idx = (t0-W):(t0+W);
        phSnip(e,:,:) = allPhase(ii, idx, :);  % [1 x nOff x nF]
        powSnip(e,:,:)= allPow(  ii, idx, :);
        envSnip(e,:)  = allGam(ii, idx);

        curTrial = chanDat.trial.data(ii,1:5:end);
        rawSnip(e,:)  = curTrial(idx);
    end

    % ---- 4) ITPC across events (BREATH-WEIGHTED) ----
    % Compute per-breath mean vector, then mean across breaths equally.
    uB = unique(evBreath(:));
    Vb = nan(numel(uB), nOff, nF);

    for bb = 1:numel(uB)
        b = uB(bb);
        m = (evBreath == b);
        if any(m)
            Vb(bb,:,:) = squeeze(mean(exp(1j*phSnip(m,:,:)), 1, 'omitnan'));
        end
    end
    itpc = abs(squeeze(mean(Vb, 1, 'omitnan'))); % [nOff x nF]

    powM = squeeze(mean(powSnip, 1, 'omitnan'));
end

% ---- save event/snippet variables ----
gammaEnv.phSnip   = phSnip;
gammaEnv.powSnip  = powSnip;
gammaEnv.rawSnip  = rawSnip;   % helpful for later ERP-style checks
gammaEnv.envSnip  = envSnip;   % gamma envelope snippet itself

gammaEnv.thr      = thr;
gammaEnv.evBreath = evBreath;
gammaEnv.evT0     = evT0;
gammaEnv.offs     = offs;

if exist('itpc','var'), gammaEnv.itpc = itpc; end
if exist('powM','var'), gammaEnv.powM = powM; end

% Also helpful for later aggregation / QC
gammaEnv.nEvents = nE;
gammaEnv.nEventsPerBreath = accumarray(evBreath(:), 1, [nBreaths 1], @sum, 0);

% --- evShoulders: nE x 5 x nF envelope values at [-1, -0.5, 0, +0.5, +1] cycles ---
shoulderCycles = [-1 -0.5 0 0.5 1];                 % 1 x 5
fs_env = S.tfInfo.fs;                               % TF sampling rate for allGam
nE = numel(evT0);
nF = numel(ampFrex);
nT = size(allGam,2);

% offsets in TF samples: 5 x nF
evShoulderSampOffsets = round( shoulderCycles(:) .* (fs_env ./ ampFrex(:).') );

% allocate
evShoulders = nan(nE, numel(shoulderCycles), nF, 'single');

for f = 1:nF
    for k = 1:numel(shoulderCycles)
        idx = evT0(:) + evShoulderSampOffsets(k,f);
        valid = (idx >= 1) & (idx <= nT);

        tmp = nan(nE,1);
        lin = sub2ind([nBreaths nT], evBreath(valid), idx(valid));
        tmp(valid) = allGam(lin);

        evShoulders(:,k,f) = single(tmp);
    end
end

gammaEnv.evShoulders = evShoulders;
gammaEnv.evShoulderCycles = shoulderCycles;
gammaEnv.evShoulderSampOffsets = evShoulderSampOffsets;
gammaEnv.evShoulderSecOffsets = shoulderCycles(:) .* (1 ./ ampFrex(:).');  % 5 x nF




% ============================================================
% Surrogate ITPC distribution (BREATH-WEIGHTED) with resampling
% from the original allPhase breath×time×freq array
% ============================================================
rng(0);                      % reproducible
gammaEnv.rngSeed = 0;

% ---- quick guards ----
nBreaths = size(allPhase,1);
nT       = size(allPhase,2);
nF       = size(allPhase,3);
nE       = numel(evT0);
nOff     = 2*W + 1;
offs     = -W:W;

if nE == 0
    error('No events (evT0) available. Lower threshold or verify detection.');
end

% ---- precompute time-index matrix for the window around each event ----
idxMat = evT0(:) + offs;  % [nE x nOff]

% ---- helper: compute ITPC given a breath vector bVec (length nE) ----
compute_itpc = @(bVec) local_compute_itpc_from_allPhase_breathWeighted(allPhase, bVec, idxMat);

% ---- 1) observed ITPC (BREATH-WEIGHTED) using original event breath assignments ----
itpc_obs = compute_itpc(evBreath(:));     % [nOff x nF]

% ---- 2) shuffles ----
itpc_shuf = nan(nOff, nF, nShuf, 'single');

for ss = 1:nShuf
    if mod(ss,100)==0, fprintf('Shuffle %d / %d\n', ss, nShuf); end

    % shuffle breath assignments across events, keep event timepoints fixed
    % null snippets are pulled from allPhase(bShuf(e), idxMat(e,:), :)
    bShuf = evBreath(randperm(nE));

    itpc_shuf(:,:,ss) = single(compute_itpc(bShuf));
end

mu  = mean(itpc_shuf, 3, 'omitnan');
sig = std(itpc_shuf, 0, 3, 'omitnan');
sig(sig < 1e-6) = 1e-6;

z_itpc = (itpc_obs - double(mu)) ./ double(sig);   % [nOff x nF]

% ---- stash results ----
gammaEnv.itpc_obs = itpc_obs;
gammaEnv.itpc_mu  = mu;
gammaEnv.itpc_sig = sig;
gammaEnv.itpc_z   = z_itpc;
gammaEnv.tim      = S.tfInfo.tim;
gammaEnv.idxMat   = idxMat;   % helpful for recomputing subsets/nulls later



% --- evShNull: nShuf x nE x 5 x nF shoulders under breath-shuffled null ---
% Recreate the same breath permutations used above by resetting RNG to the same seed.
shoulderCycles = gammaEnv.evShoulderCycles;
evShoulderSampOffsets = gammaEnv.evShoulderSampOffsets;

evShNull = nan(nShuf, nE, numel(shoulderCycles), nF, 'single');

rng(gammaEnv.rngSeed);

for ss = 1:nShuf
    bShuf = evBreath(randperm(nE));  % shuffle breaths across events (event times fixed)

    for f = 1:nF
        for k = 1:numel(shoulderCycles)
            idx = evT0(:) + evShoulderSampOffsets(k,f);
            valid = (idx >= 1) & (idx <= nT);

            tmp = nan(nE,1);
            lin = sub2ind([nBreaths nT], bShuf(valid), idx(valid));
            tmp(valid) = allGam(lin);

            evShNull(ss,:,k,f) = single(tmp(:)).';
        end
    end
end

gammaEnv.evShNull = evShNull;


muNull  = squeeze(mean(gammaEnv.evShNull, 1, 'omitnan'));        % [nE x 5 x nF]
sigNull = squeeze(std(gammaEnv.evShNull, 0, 1, 'omitnan'));      % [nE x 5 x nF]
sigNull(sigNull < 1e-6) = 1e-6;

gammaEnv.evShoulders_z = (double(gammaEnv.evShoulders) - double(muNull)) ./ double(sigNull); % [nE x 5 x nF]


% Bootstrap ERP with at most 1 event per breath, and 30 events per bootstrap
% Requires:
%   gammaEnv.envSnip   : [nE x nT]
%   gammaEnv.evBreath  : [nE x 1] breath index for each event

envSnip = gammaEnv.envSnip;
evBreath = gammaEnv.evBreath(:);

nBoot = 1000;
nPerBoot = 30;

[nE, nT] = size(envSnip);

% --- choose one representative event per breath (randomly if multiple) ---
uB = unique(evBreath);
repIdx = nan(numel(uB),1);

for bb = 1:numel(uB)
    b = uB(bb);
    m = find(evBreath == b);
    repIdx(bb) = m(randi(numel(m)));   % pick 1 event for this breath
end

envRep = envSnip(repIdx, :);           % [nBreathsWithEvents x nT]
nAvail = size(envRep, 1);

if nAvail < nPerBoot
    warning('Only %d breaths available; bootstraps will sample with replacement to reach %d.', nAvail, nPerBoot);
end

% --- bootstrap ---
erp_boot_all = nan(nBoot, nT);

for b = 1:nBoot
    idx = randi(nAvail, [nPerBoot 1]);                     % 30 events (breaths) with replacement
    erp_boot_all(b,:) = mean(envRep(idx,:), 1, 'omitnan'); % mean ERP for this bootstrap
end

erp_boot_mean = median(erp_boot_all, 1, 'omitnan');
erp_boot_se   = std(erp_boot_all, 0, 1, 'omitnan');
erp_boot_ci   = prctile(erp_boot_all, [2.5 97.5], 1);

% stash if desired
gammaEnv.erp_boot_mean = erp_boot_mean;
gammaEnv.erp_boot_se   = erp_boot_se;
gammaEnv.erp_boot_ci   = erp_boot_ci;
% gammaEnv.erp_boot_all  = erp_boot_all;  % optional (large)



end

% =====================================================================
% Local function: BREATH-WEIGHTED ITPC computed by resampling snippets
% from allPhase(breath, time, freq) using the provided bVec + idxMat.
% =====================================================================
function itpc = local_compute_itpc_from_allPhase_breathWeighted(allPhase, bVec, idxMat)
% allPhase: [nBreaths x nT x nF]
% bVec:     [nE x 1] breath index per event
% idxMat:   [nE x nOff] time indices per event window

nF   = size(allPhase,3);
nOff = size(idxMat,2);

uB = unique(bVec(:));
Vb = nan(numel(uB), nOff, nF);

for bb = 1:numel(uB)
    b = uB(bb);
    eIdx = find(bVec == b);
    if isempty(eIdx), continue; end

    idx = idxMat(eIdx,:); % [k x nOff]

    for f = 1:nF
        ph2 = allPhase(b,:,f);          % [1 x nT]
        ph  = ph2(idx);                 % [k x nOff]
        Vb(bb,:,f) = mean(exp(1j*ph), 1, 'omitnan'); % [1 x nOff]
    end
end

itpc = abs(squeeze(mean(Vb, 1, 'omitnan'))); % [nOff x nF]
end