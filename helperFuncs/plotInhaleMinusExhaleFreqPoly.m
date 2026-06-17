function [hFig, subInfo, sortedBreathIdxCell, diffCell] = ...
    plotInhaleMinusExhaleFreqPoly(idxVec, subjectInfo, dataCell, ...
    bgCol, labCol, dotColors, binWidth)
% plotInhaleMinusExhaleFreqPoly
%
% For each selected cell in dataCell:
%   - assumes data are [breath x block(5) x frequency]
%   - for each of the 5 respiratory epochs, finds:
%         (1) the row-wise peak power in 25-60 Hz
%         (2) the frequency at which that peak power occurs
%   - defines final inhale / exhale values as:
%         inhale = frequency from whichever of epochs 1:2 has larger peak power
%         exhale = frequency from whichever of epochs 3:4 has larger peak power
%   - computes a per-breath difference score:
%         diffHz = inhaleHz - exhaleHz
%   - plots a ggplot-style frequency polygon for each participant
%   - returns, for each participant instance, the ORIGINAL breath indices ranked from
%     largest frequency shift to smallest
%
% IMPORTANT:
%   "largest frequency shift" is interpreted here as largest ABSOLUTE
%   inhale-exhale frequency difference, i.e. sorted by abs(diffHz) descending.
%
% PLOTTING CHANGE:
%   If multiple retained rows belong to the same subject name, their breath-level
%   diffHz values are concatenated into a single plotted distribution.
%   Outputs remain split apart per retained row / participant instance.
%
% Y-AXIS CHANGE:
%   Each plotted subject distribution is shown as proportion within bin:
%       counts / total_counts_for_that_subject
%   so the proportions across bins sum to 1 for each plotted subject.
%
% Outputs:
%   hFig               : figure handle
%   subInfo            : summary table, one row per retained participant instance
%   sortedBreathIdxCell: cell array, one cell per retained participant instance;
%                        each cell contains original breath indices sorted
%                        from largest |inhale-exhale| shift to smallest
%   diffCell           : cell array, one cell per retained participant instance;
%                        each cell contains the signed inhale-exhale
%                        difference values for the retained breaths, in the
%                        same order as sortedBreathIdxCell
%
% Assumption:
%   inhale uses epochs 1:2
%   exhale uses epochs 3:4
%   pause (epoch 5) is calculated but not used
%
% Optional:
%   dotColors : [nSelected x 3] RGB colors, one per selected cell
%   binWidth  : histogram bin width in Hz (default = 1)

frex = logspace(log10(.1),log10(200),300);

% --- pull selected cells + names ---
dup1Cells = dataCell(idxVec);
subNames  = string(subjectInfo(idxVec,8));
subNames  = strtrim(subNames);

nCell = numel(dup1Cells);

% optional colors
if nargin < 6 || isempty(dotColors)
    dotColors = lines(nCell);
elseif size(dotColors,1) < nCell
    error('dotColors must have at least as many rows as selected cells.');
end

% optional bin width
if nargin < 7 || isempty(binWidth)
    binWidth = 1;
end

% peak frequency band
fMask = frex >= 25 & frex <= 60;
fBand = frex(fMask);

% summary table
subInfo = table( ...
    (1:nCell)', ...
    subNames(:), ...
    nan(nCell,1), ...
    nan(nCell,1), ...
    nan(nCell,1), ...
    nan(nCell,1), ...
    'VariableNames', {'cellIdx','subName','nBreaths','meanDiffHz','sdDiffHz','medianDiffHz'});

% first pass: compute per-instance difference scores
diffCellAll = cell(nCell,1);
sortedBreathIdxCellAll = cell(nCell,1);
allDiffs = [];

for i = 1:nCell
    tmp = dup1Cells{i};

    % expect [breath x 5 x frequency]
    if ndims(tmp) ~= 3
        warning('Skipping cell %d: expected 3D array [breath x block x freq].', i);
        continue
    end
    if size(tmp,2) < 5
        warning('Skipping cell %d: second dimension has fewer than 5 levels.', i);
        continue
    end

    nBreathOrig = size(tmp,1);
    origBreathIdx = (1:nBreathOrig)';

    % keep only first 5 levels of dim 2, and 25-60 Hz in dim 3
    tmp = tmp(:,1:5,fMask);   % [breath x 5 x nBand]

    peakHz  = nan(nBreathOrig, 5);
    peakPow = nan(nBreathOrig, 5);

    % row-wise peak power and corresponding frequency for each epoch
    for b = 1:5
        cur = squeeze(tmp(:,b,:));   % [breath x nBand]
        if isvector(cur)
            cur = cur(:)';
        end

        goodRow = any(isfinite(cur), 2);
        if any(goodRow)
            [pkPow, pkIdx] = max(cur(goodRow,:), [], 2);
            peakPow(goodRow,b) = pkPow;
            peakHz(goodRow,b)  = fBand(pkIdx);
        end
    end

    % keep only breaths with all 5 epoch peaks defined
    keepRow = all(isfinite(peakHz), 2) & all(isfinite(peakPow), 2);
    peakHz  = peakHz(keepRow,:);
    peakPow = peakPow(keepRow,:);
    keepBreathIdx = origBreathIdx(keepRow);

    if isempty(peakHz)
        continue
    end

    % inhale: choose frequency from epochs 1:2 with larger peak power
    [~, inhSel] = max(peakPow(:,1:2), [], 2);   % 1 or 2
    inhaleHz = peakHz(sub2ind(size(peakHz), (1:size(peakHz,1))', inhSel));

    % exhale: choose frequency from epochs 3:4 with larger peak power
    [~, exSelLocal] = max(peakPow(:,3:4), [], 2); % 1 or 2
    exSel = exSelLocal + 2;                       % 3 or 4
    exhaleHz = peakHz(sub2ind(size(peakHz), (1:size(peakHz,1))', exSel));

    diffHz = inhaleHz - exhaleHz;
    goodDiff = isfinite(diffHz);

    diffHz = diffHz(goodDiff);
    keepBreathIdx = keepBreathIdx(goodDiff);

    if isempty(diffHz)
        continue
    end

    % sort ORIGINAL breath indices by largest absolute shift to smallest
    [~, sortOrd] = sort(abs(diffHz), 'descend');
    sortedBreathIdxCellAll{i} = keepBreathIdx(sortOrd);
    diffCellAll{i} = diffHz(sortOrd);

    subInfo.nBreaths(i)     = numel(diffHz);
    subInfo.meanDiffHz(i)   = mean(diffHz, 'omitnan');
    subInfo.sdDiffHz(i)     = std(diffHz, 'omitnan');
    subInfo.medianDiffHz(i) = median(diffHz, 'omitnan');

    allDiffs = [allDiffs; diffHz(:)];
end

% drop empty rows
keepSub = subInfo.nBreaths > 0;
subInfo = subInfo(keepSub,:);
diffCell = diffCellAll(keepSub);
sortedBreathIdxCell = sortedBreathIdxCellAll(keepSub);
subColors = dotColors(keepSub,:);

% --- combine identical subjects for plotting only ---
plotSubNames = unique(subInfo.subName, 'stable');
nPlotSub = numel(plotSubNames);

plotDiffCell = cell(nPlotSub,1);
plotColors   = nan(nPlotSub,3);

for p = 1:nPlotSub
    thisRows = find(subInfo.subName == plotSubNames(p));

    tmpDiff = cell(numel(thisRows),1);
    for r = 1:numel(thisRows)
        tmpDiff{r} = diffCell{thisRows(r)}(:);
    end
    plotDiffCell{p} = vertcat(tmpDiff{:});

    % use first occurrence color for plotting
    plotColors(p,:) = subColors(thisRows(1),:);
end

% --- figure ---
hFig = figure('Color', bgCol, 'Position', [50 50 770*.9 490*.9], 'InvertHardcopy', 'off');
ax = axes('Parent', hFig);
hold(ax, 'on');

ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.YColor     = labCol;
ax.TickDir    = 'out';
ax.TickLength = [0.018 0.018];
box(ax, 'off');

if isempty(allDiffs)
    ylabel('Proportion within bin', ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', labCol);

    xlabel('Inhale - exhale peak frequency difference (Hz)', ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', labCol);

    % title('Distribution of inhale - exhale peak frequency differences', ...
    %     'FontSize', 20, ...
    %     'FontWeight', 'bold', ...
    %     'FontName', 'Dotum', ...
    %     'Color', labCol);

    xlim([-5 5])
    ylim([0 1])
    set(ax, 'Layer', 'top')
    return
end

% common bins across all participant instances / breaths
xMin = floor(min(allDiffs) / binWidth) * binWidth;
xMax = ceil(max(allDiffs) / binWidth) * binWidth;
if xMin == xMax
    xMin = xMin - binWidth;
    xMax = xMax + binWidth;
end
edges = xMin:binWidth:xMax;
if numel(edges) < 2
    edges = [xMin, xMin + binWidth];
end
centers = edges(1:end-1) + diff(edges(1:2))/2;

maxProp = 0;

% plot frequency polygon for each UNIQUE subject, after combining rows
for i = 1:nPlotSub
    thisDiff = plotDiffCell{i};
    thisCol  = plotColors(i,:);

    counts = histcounts(thisDiff, edges);

    if sum(counts) > 0
        props = counts ./ sum(counts);
    else
        props = zeros(size(counts));
    end

    maxProp = max(maxProp, max(props));

    plot(ax, centers, props, '-', ...
        'Color', thisCol, ...
        'LineWidth', 10, ...
        'HandleVisibility', 'off');
end

% zero reference line
plot(ax, [0 0], [0 max(maxProp,1e-6)], '--', ...
    'Color', labCol, ...
    'LineWidth', 7, ...
    'HandleVisibility', 'off');

ylabel('Proportion within bin', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

xlabel('Inhale - exhale peak frequency difference (Hz)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

% title('Distribution of inhale - exhale peak frequency differences', ...
%     'FontSize', 20, ...
%     'FontWeight', 'bold', ...
%     'FontName', 'Dotum', ...
%     'Color', labCol);

set(ax, 'Layer', 'top')

xlim([-7.5 15])
ylim([0 max(0.05, maxProp * 1.08)])

end