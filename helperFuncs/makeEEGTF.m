function fig = makeEEGTF(chanDat, taskVec, conds, opts)
%MAKEEEGTF  Plot mean time-frequency (powZ) heatmaps by condition.
%
%   fig = makeEEGTF(chanDat, taskVec, conds)
%   fig = makeEEGTF(chanDat, taskVec, conds, opts)
%
% INPUTS
%   chanDat.tf.powZ : [nBreath x nTime x nFreq] z-scored power TF
%   chanDat.tf.tVec : [1 x nTime] time vector (optional)
%   chanDat.tf.frex : [1 x nFreq] frequency vector (optional)
%   chanDat.use     : [nBreath x 1] logical/0-1
%   taskVec         : string array length nBreath with condition labels
%   conds           : string array, e.g. ["audio","focus","shadow"]
%   opts (optional struct):
%       .clim       (default [-5 5])
%       .colormap   (default parula)
%       .minN       (default 5) minimum breaths to plot a condition
%
% OUTPUT
%   fig : figure handle

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'clim')     || isempty(opts.clim),     opts.clim = [-5 5]; end
if ~isfield(opts,'colormap') || isempty(opts.colormap), opts.colormap = parula; end
if ~isfield(opts,'minN')     || isempty(opts.minN),     opts.minN = 5; end

powZ = chanDat.tf.powZ;  % [breath x time x freq]
assert(ndims(powZ)==3, 'chanDat.tf.powZ must be [breaths x time x frequency].');

[nBreath, nTime, nFreq] = size(powZ);

useVec = logical(chanDat.use(:));
assert(numel(useVec)==nBreath, 'chanDat.use length must match size(powZ,1).');

taskVec = string(taskVec(:));
assert(numel(taskVec)==nBreath, 'taskVec length must match size(powZ,1).');

conds = string(conds(:))';


tVec = 1:nTime;


fVec = double(chanDat.tf.frex(:))';


fig = figure('Color','w', 'visible', false, 'position', [0,0,1200, 700]);
tlo = tiledlayout(fig, 2, numel(conds), 'Padding','compact', 'TileSpacing','compact');

for c = 1:numel(conds)
    ax = nexttile(tlo, c);
    m = useVec & (taskVec == conds(c));

    if nnz(m) < opts.minN
        imagesc(ax, [],[], nan(numel(fVec), numel(tVec)));
        axis(ax,'xy');
        
        title(ax, sprintf('%s (n=%d)', conds(c), nnz(m)));
        xlabel(ax,'Time'); ylabel(ax,'Freq');
        caxis(ax, opts.clim);
        colormap(ax, opts.colormap);
        colorbar(ax);
        continue
    end

    fb = mean(1 ./ chanDat.behDat.length(m));

    [~, fbi] = min(abs(fb - fVec)); 
    tfMean = squeeze(mean(powZ(m,:,:), 1, 'omitnan')); % [nTime x nFreq]
    tfMean = tfMean.';                                 % [nFreq x nTime]

    imagesc(ax, [], [], tfMean);
    axis(ax,'xy');
    yline(fbi)
    yticks(20:20:150)
    title(ax, sprintf('%s (n=%d)', conds(c), nnz(m)));
    yticklabels(round(fVec(20:20:150)))
    xlabel(ax,'Time'); ylabel(ax,'Freq');
    caxis(ax, opts.clim);
    colormap(ax, opts.colormap);

    cb = colorbar(ax);
    cb.Label.String = 'powZ';
end


% --- Colors to match earlier scatter/condition scheme ---
cols = lines(numel(conds));

freqBreaks = [2, 8; ...
              8, 14;...
              14,30];
freqLabs = {'theta', 'alpha', 'beta'};

for fi = 1:3
    ax = nexttile(tlo, fi+3);
    cla(ax); hold(ax,'on');

    bandMask = (fVec >= freqBreaks(fi,1)) & (fVec < freqBreaks(fi,2));

    % Robust median across time (10:20) and freq band for each breath
    A = chanDat.tf.powZ(:,:, bandMask);  % [nBreath x nTimeWin x nFband]
    medVals = squeeze(max(max(A, [],2, 'omitnan'), [], 3, 'omitnan')); % [nBreath x 1]
    medVals = double(medVals(:));

    % Bin edges shared across conditions for this band
    mAll = useVec(:) & ~ismissing(taskVec) & taskVec ~= "" & isfinite(medVals);
    valsAll = medVals(mAll);

    if isempty(valsAll)
        title(ax, sprintf('%s power (no data)', freqLabs{fi}));
        box(ax,'off'); hold(ax,'off');
        continue
    end

    lo = prctile(valsAll, 1);
    hi = prctile(valsAll, 95);
    if lo == hi
        lo = lo - 1; hi = hi + 1;
    end
    nBins = 10;
    binEdges = linspace(lo, hi, nBins+1);

    % Plot histograms per condition
    for c = 1:numel(conds)
        m = useVec(:) & (taskVec == conds(c)) & isfinite(medVals);

        if nnz(m) < opts.minN
            continue
        end

        histogram(ax, medVals(m), binEdges, ...
            'DisplayStyle','stairs', ...
            'EdgeColor', cols(c,:), ...
            'LineWidth', 2, ...
            'HandleVisibility','off', ...
            'normalization', 'probability');
    end

    title(ax, sprintf('%s band (median powZ)', freqLabs{fi}));
    xlabel(ax, 'median powZ (time 10:20, band)');
    ylabel(ax, 'Count');
    box(ax,'off');
    grid(ax,'on');
    hold(ax,'off');
end














% ============================================================
% Layout: 6 x 3 tiledlayout
%   Rows 1-3: heatmaps (audio/focus/shadow), each spans all 3 columns
%   Rows 4-6: histograms by condition (audio/focus/shadow)
%             Col 1=theta, Col 2=alpha, Col 3=beta
%             Same x-lims (and optional y-lims) within each band across rows
% ============================================================

fig = figure('Color','w', 'visible', false, 'position', [0,0,1200, 900]);
tlo = tiledlayout(fig, 6, 3, 'Padding','compact', 'TileSpacing','compact');

cols = lines(numel(conds));  % condition colors (audio/focus/shadow)

% -------------------- TOP: heatmaps (UNCHANGED logic; new tile placement) --------------------
for c = 1:numel(conds)
    ax = nexttile(tlo, c, [3 1]);
    m  = useVec & (taskVec == conds(c));
    freqMask = fVec>=2 & fVec<=30; 
    if nnz(m) < opts.minN
        imagesc(ax, [], [], nan(sum(freqMask), numel(tVec)));
        axis(ax,'xy');
        title(ax, sprintf('%s (n=%d)', conds(c), nnz(m)));
        xlabel(ax,'Time'); ylabel(ax,'Freq');
        caxis(ax, opts.clim);
        colormap(ax, opts.colormap);
        colorbar(ax);
        continue
    end

    fb = mean(1 ./ chanDat.behDat.length(m));
    [~, fbi] = min(abs(fb - fVec));

    tfMean = squeeze(mean(powZ(m,:,freqMask), 1, 'omitnan')); % [nTime x nFreq]
    tfMean = tfMean.';                                 % [nFreq x nTime]

    imagesc(ax, [], [], tfMean);
    axis(ax,'xy');
    % yline(ax, fbi);
    
    yticks(ax, 10:10:60);
    curFreq = fVec(freqMask); 
    yticklabels(ax, round(curFreq(10:10:60)));

    bands = [8,14];
    [~, fi] = min(abs(bands(1) - curFreq));
    yline(fi)
    [~, fi] = min(abs(bands(2) - curFreq));
    yline(fi)

    title(ax, sprintf('%s (n=%d)', conds(c), nnz(m)));
    xlabel(ax,'Time'); ylabel(ax,'Freq');

    caxis(ax, opts.clim);
    colormap(ax, opts.colormap);

    cb = colorbar(ax);
    cb.Label.String = 'powZ';
end

% -------------------- BOTTOM: histograms (separate axes per condition) --------------------
freqBreaks = [2, 8; 8, 14; 14, 30];
freqLabs   = {'theta','alpha','beta'};

nB = size(chanDat.tf.powZ, 1);
bandVals = nan(nB, 3);

% breath-wise band summary (edit this summary if you want mean/max instead of median)
for fi = 1:3
    bandMask = (fVec >= freqBreaks(fi,1)) & (fVec < freqBreaks(fi,2));
    A = chanDat.tf.powZ(:, :, bandMask);  % [nBreath x timeWin x fBand]
    bandVals(:,fi) = squeeze(median(squeeze(prctile(A, 90, 2)),2,'omitnan')); % [nBreath x 1]
end

% shared x-lims + bins per band (across all conditions)
xlimBand = nan(3,2);
binEdges = cell(1,3);
nBins = 12;

mBase = useVec(:) & ~ismissing(taskVec) & taskVec ~= "";

for fi = 1:3
    v = bandVals(:,fi);
    vv = v(mBase & isfinite(v));
    if isempty(vv)
        xlimBand(fi,:) = [0 1];
    else
        lo = prctile(vv, 1);
        hi = prctile(vv, 95);
        if lo == hi, lo = lo-1; hi = hi+1; end
        xlimBand(fi,:) = [lo hi];
    end
    binEdges{fi} = linspace(xlimBand(fi,1), xlimBand(fi,2), nBins+1);
end

% optional: shared y-lims per band when using probability normalization
yMaxBand = zeros(1,3);
for fi = 1:3
    v = bandVals(:,fi);
    for c = 1:numel(conds)
        mc = mBase & (taskVec == conds(c)) & isfinite(v);
        if nnz(mc) < opts.minN, continue; end
        p = histcounts(v(mc), binEdges{fi}, 'Normalization','probability');
        yMaxBand(fi) = max(yMaxBand(fi), max(p));
    end
end

% plot: rows 4-6 are conditions; cols 1-3 are theta/alpha/beta
for c = 1:numel(conds)
    for fi = 1:3
        ax = nexttile(tlo, (3 + (c-1))*3 + fi);

        v  = bandVals(:,fi);
        mc = mBase & (taskVec == conds(c)) & isfinite(v);

        cla(ax); hold(ax,'on');

        if nnz(mc) < opts.minN
            title(ax, sprintf('%s: %s (n=%d)', conds(c), freqLabs{fi}, nnz(mc)));
            xlim(ax, xlimBand(fi,:));
            box(ax,'off');
            hold(ax,'off');
            continue
        end

        histogram(ax, v(mc), binEdges{fi}, ...
            'Normalization','probability', ...
            'FaceColor', cols(c,:), ...
            'EdgeColor','none', ...
            'FaceAlpha', 0.85);

        xlim(ax, xlimBand(fi,:));
        if yMaxBand(fi) > 0
            ylim(ax, [0 yMaxBand(fi)*1.05]);
        end

        title(ax, sprintf('%s: %s (n=%d)', conds(c), freqLabs{fi}, nnz(mc)));

        % tidy labeling (reduce clutter)
        if fi == 1
            ylabel(ax, 'Probability');
        else
            ax.YTickLabel = [];
        end
        if c == numel(conds)
            xlabel(ax, 'band powZ');
        else
            ax.XTickLabel = [];
        end

        grid(ax,'on');
        box(ax,'off');
        hold(ax,'off');
    end
end


end