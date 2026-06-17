function [hFig, subInfo] = plotInhaleMinusExhalePowerFreqPoly(idxVec, subjectInfo, dataCell, ...
    bgCol, labCol, dotColors, binWidth)
% plotInhaleMinusExhalePowerFreqPoly
%
% For each selected cell in dataCell:
%   - assumes data are [breath x block(5) x frequency]
%   - for each of the 5 respiratory epochs, finds the row-wise peak power
%     in the 25-60 Hz band
%   - defines final inhale / exhale peak-power values as:
%         inhale = larger peak power from epochs 1:2
%         exhale = larger peak power from epochs 3:4
%   - computes a per-breath difference score:
%         diffPow = inhalePow - exhalePow
%   - plots a ggplot-style frequency polygon for each participant in their color
%
% Styling is matched to the other plotting functions.
%
% Assumption:
%   inhale uses epochs 1:2
%   exhale uses epochs 3:4
%   pause (epoch 5) is calculated but not used
%
% Optional:
%   dotColors : [nSelected x 3] RGB colors, one per selected cell
%   binWidth  : histogram bin width in power units (default = 0.25)

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
    binWidth = 0.25;
end

% peak frequency band
fMask = frex >= 25 & frex <= 60;

% summary table
subInfo = table( ...
    (1:nCell)', ...
    subNames(:), ...
    nan(nCell,1), ...
    nan(nCell,1), ...
    nan(nCell,1), ...
    nan(nCell,1), ...
    'VariableNames', {'cellIdx','subName','nBreaths','meanDiffPow','sdDiffPow','medianDiffPow'});

% first pass: compute per-subject difference scores
diffCell = cell(nCell,1);
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

    % keep only first 5 levels of dim 2, and 25-60 Hz in dim 3
    tmp = tmp(:,1:5,fMask);   % [breath x 5 x nBand]

    nBreath = size(tmp,1);
    peakPow = nan(nBreath, 5);

    % row-wise peak power for each epoch
    for b = 1:5
        cur = squeeze(tmp(:,b,:));   % [breath x nBand]
        if isvector(cur)
            cur = cur(:)';
        end

        goodRow = any(isfinite(cur), 2);
        if any(goodRow)
            pkPow = max(cur(goodRow,:), [], 2);
            peakPow(goodRow,b) = pkPow;
        end
    end

    % keep only breaths with all 5 epoch peaks defined
    keepRow = all(isfinite(peakPow), 2);
    peakPow = peakPow(keepRow,:);

    if isempty(peakPow)
        continue
    end

    % inhale / exhale peak-power values
    inhalePow = max(peakPow(:,1:2), [], 2);
    exhalePow = max(peakPow(:,3:4), [], 2);

    diffPow = inhalePow - exhalePow;
    diffPow = diffPow(isfinite(diffPow));

    if isempty(diffPow)
        continue
    end

    diffCell{i} = diffPow;

    subInfo.nBreaths(i)      = numel(diffPow);
    subInfo.meanDiffPow(i)   = mean(diffPow, 'omitnan');
    subInfo.sdDiffPow(i)     = std(diffPow, 'omitnan');
    subInfo.medianDiffPow(i) = median(diffPow, 'omitnan');

    allDiffs = [allDiffs; diffPow(:)];
end

% drop empty rows
keepSub = subInfo.nBreaths > 0;
subInfo = subInfo(keepSub,:);
diffCell = diffCell(keepSub);
subColors = dotColors(keepSub,:);

% --- figure ---
hFig = figure('Color', bgCol, 'Position', [50 50 1100 700], 'InvertHardcopy', 'off');
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
    ylabel('Count', ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', labCol);

    xlabel('Inhale - exhale peak power difference', ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', labCol);

    title('Distribution of inhale - exhale peak power differences', ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', labCol);

    xlim([-1 1])
    ylim([0 1])
    set(ax, 'Layer', 'top')
    return
end

% common bins across all participants
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

maxCount = 0;

% plot frequency polygon for each subject
for i = 1:numel(diffCell)
    thisDiff = diffCell{i};
    thisCol  = subColors(i,:);

    counts = histcounts(thisDiff, edges);
    maxCount = max(maxCount, max(counts));

    plot(ax, centers, counts, '-', ...
        'Color', thisCol, ...
        'LineWidth', 2.4, ...
        'HandleVisibility', 'off');
end

% zero reference line
plot(ax, [0 0], [0 max(maxCount,1)], '--', ...
    'Color', labCol, ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off');

ylabel('Count', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

xlabel('Inhale - exhale peak power difference', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

title('Distribution of inhale - exhale peak power differences', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

set(ax, 'Layer', 'top')

xlim([-7.5 15])
ylim([0 max(1, ceil(maxCount * 1.08))])

end