function allChannel = makeAllChannel()
% makeAllChannel
%   Groups chanDat files by stem (filename without final _NN),
%   picks the best-gamma channel per stem using the same logic as your pipeline:
%       score = robust_z(prevPooled) + robust_z(promPooled)
%   and returns a struct array with selected fields from each winning chanDat.
%
% Usage:
%   allChannel = makeAllChannel();

% ---------- USER PATH (matches your pipeline) ----------
procDir = "R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\CHANDAT_processed\";

% ---------- Selection params (matches your pipeline) ----------
gammaBandHz          = [25 60]; %#ok<NASGU>  % included for signature parity; not used in prom calc (matches your code)
baselineBandHz       = [25 58];
excludeHzAroundPeak  = 5;

% ---------- Find .mat files (top level only) ----------
files = dir(fullfile(procDir, "*.mat"));
files = files(~[files.isdir]);
test = cellfun(@(x) sum(strfind(x, 'EEG'))>0, {files.name}); 
files(test) = []; 
if isempty(files)
    error("No .mat files found in %s", procDir);
end



% ---------- Group key = filename minus final _NN ----------
keys = strings(numel(files),1);
for i = 1:numel(files)
    keys(i) = stripChanSuffix(files(i).name);
end
[uKeys, ~, gidx] = unique(keys, 'stable');
nGroups = numel(uKeys);

% ---------- Preallocate output ----------
tmpl = struct( ...
    'task', [], 'sessID', [], 'subID', [], 'sessNum', [], 'type', [], ...
    'use', [], 'reasonEliminate', [], 'tim', [], 'tf', [], 'behDat', [], 'fooof', [], ...
    ... % helpful provenance (extra; remove if you truly want ONLY the requested fields)
    'groupKey', [], 'winnerFilePath', [], 'winnerFileName', [], 'chi', [], 'chanLabel', [] ...
);
allChannel = repmat(tmpl, nGroups, 1);

fprintf("Found %d stem-groups.\n", nGroups);

% ---------- Loop over stem-groups ----------
for gi = 1:nGroups
    gi
    idx = find(gidx == gi);
    groupFiles = files(idx);
    nC = numel(groupFiles);

    prevPooled = nan(nC,1);
    promPooled = nan(nC,1);

    % summarize each channel
    for ci = 1:nC
        fpath = fullfile(groupFiles(ci).folder, groupFiles(ci).name);
        S = load(fpath, "chanDat");
        cd = S.chanDat;

        useVec = cd.use(:) == 1;

        % prevalence pooled across epochs
        G = double(cd.fooof.gamma_peaks); % breaths x epochs (NaN = no peak)
        nEpochs = size(G,2);
        den = sum(useVec);

        prev_by_epoch = nan(1,nEpochs);
        for e = 1:nEpochs
            prev_by_epoch(e) = sum(useVec & isfinite(G(:,e))) / max(den,1);
        end
        prevPooled(ci) = mean(prev_by_epoch, 'omitnan');

        % prominence pooled across breaths x epochs (your existing definition)
        promPooled(ci) = gamma_prom_pooled(cd, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak);
    end

    % winner selection (same robust z logic; no plotting)
    zPrev  = robust_z(prevPooled);
    zProm  = robust_z(promPooled);
    score  = zPrev + zProm;

    [~, ord] = sort(score, 'descend', 'MissingPlacement', 'last');
    win = ord(1);

    winnerPath = fullfile(groupFiles(win).folder, groupFiles(win).name);

    % load winner and extract requested fields
    Sw = load(winnerPath, "chanDat");
    cw = Sw.chanDat;

    allChannel(gi).task            = cw.task;
    allChannel(gi).sessID          = cw.sessID;
    allChannel(gi).subID           = cw.subID;
    allChannel(gi).sessNum         = cw.sessNum;
    allChannel(gi).type            = cw.type;
    allChannel(gi).use             = cw.use;
    allChannel(gi).reasonEliminate = cw.reasonEliminate;
    allChannel(gi).tim             = cw.tim;
    allChannel(gi).tf              = cw.tf;
    allChannel(gi).behDat          = cw.behDat;
    allChannel(gi).fooof           = cw.fooof;

    % provenance (optional)
    allChannel(gi).groupKey        = uKeys(gi);
    allChannel(gi).winnerFilePath  = string(winnerPath);
    allChannel(gi).winnerFileName  = string(groupFiles(win).name);
    allChannel(gi).chi             = cw.chi;
    allChannel(gi).chanLabel       = chan_label(cw);

    fprintf("(%4d/%4d) %s  ->  %s (score=%.3f)\n", gi, nGroups, uKeys(gi), groupFiles(win).name, score(win));
end

end % makeAllChannel


% ========================= Local helpers =========================

function s = stripChanSuffix(fname)
% Remove final _NN.mat  (e.g., abc_def_03.mat -> abc_def)
s = regexprep(string(fname), "_\d\d\.mat$", "");
end

function m = robust_mad(x)
x = x(:);
x = x(isfinite(x));
if isempty(x), m = NaN; return; end
med = median(x);
m = median(abs(x - med));
end

function z = robust_z(x)
x = x(:);
med = median(x, 'omitnan');
m = robust_mad(x);
if ~isfinite(m) || m==0
    z = (x - med);
else
    z = (x - med) ./ (1.4826*m);
end
end

function promPooled = gamma_prom_pooled(chanDat, useVec, gammaBandHz, baselineBandHz, excludeHzAroundPeak) %#ok<INUSD>
% Matches your compute_gamma_prominence_full -> pooled median logic

Gdet = double(chanDat.fooof.gamma_peaks);        % NaN if no peak
Gmax = double(chanDat.fooof.gamma_peak_freq);    % fallback
flat = double(chanDat.fooof.spectra_flat_log10); % breaths x epochs x frex
frex = double(chanDat.tf.frex(:));

Guse = Gdet;
missing = ~isfinite(Guse);
Guse(missing) = Gmax(missing);

maskBaseBand = frex >= baselineBandHz(1) & frex <= baselineBandHz(2);

[nBreaths, nEpochs] = size(Guse);
promAll = nan(0,1);

for b = 1:nBreaths
    if ~useVec(b), continue; end
    for e = 1:nEpochs
        fpk = Guse(b,e);
        if ~isfinite(fpk), continue; end

        [~, fiPk] = min(abs(frex - fpk));

        maskEx = maskBaseBand & ~(frex >= (fpk-excludeHzAroundPeak) & frex <= (fpk+excludeHzAroundPeak));
        base = median(flat(b,e,maskEx), 3, 'omitnan');
        if ~isfinite(base), continue; end

        p = flat(b,e,fiPk) - base;
        if isfinite(p)
            promAll(end+1,1) = p; %#ok<AGROW>
        end
    end
end

promPooled = median(promAll, 'omitnan');
end

function s = chan_label(chanDat)
try
    s = string(chanDat.labels{chanDat.chi});
catch
    s = "chan" + string(chanDat.chi);
end
end