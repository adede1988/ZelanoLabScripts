function [hFig, subInfo] = plotPeakFreq_bySubjectBlocks(idxVec, subjectInfo, dataCell, ...
    bgCol, labCol, dotColors)
% plotPeakFreq_bySubjectBlocks
%
% For each selected cell in dataCell:
%   - assumes data are [breath x block(5) x frequency]
%   - first finds, for each of the 5 respiratory epochs:
%         (1) the row-wise peak power in 25-60 Hz
%         (2) the frequency at which that peak power occurs
%   - then collapses those 5 epoch-wise values into:
%         inhale = frequency from whichever of epochs 1:2 has the larger peak power
%         exhale = frequency from whichever of epochs 3:4 has the larger peak power
%   - then 0-centers each subject by subtracting that subject's overall mean
%     across all final inhale and exhale points
%   - plots jittered dots with faint within-breath linking lines
%   - overlays a thicker line for each subject mean
%
% Assumption here:
%   "two early epochs" = epochs 1:2
%   "two late epochs"  = epochs 3:4
%   pause (epoch 5) is calculated but not used in inhale/exhale selection
%
% Colors indicate which dup1Cells index the breath came from.
%
% Optional:
%   dotColors : [nSelected x 3] RGB colors, one per selected cell

frex = logspace(log10(.1),log10(200),300);

phaseNames = {'inhale', 'exhale'};
xBase = 1:2;

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
    nan(nCell,1), ...
    nan(nCell,1), ...
    'VariableNames', {'cellIdx','subName','nBreaths','subjectGrandMeanHz', ...
                      'meanInhaleHz','sdInhaleHz','meanExhaleHz','sdExhaleHz'});

% --- figure ---
hFig = figure('Color', bgCol, 'Position', [50 50 500 700], 'InvertHardcopy', 'off');
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

allPeakVals = [];

% --- loop over selected cells ---
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
    peakHz  = nan(nBreath, 5);
    peakPow = nan(nBreath, 5);

    % row-wise peak power and its frequency for each block level
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

    % keep only breaths with all 5 block peaks defined
    keepRow = all(isfinite(peakHz), 2) & all(isfinite(peakPow), 2);
    peakHz  = peakHz(keepRow,:);
    peakPow = peakPow(keepRow,:);

    if isempty(peakHz)
        continue
    end

    % inhale: choose the frequency from epochs 1:2 with larger peak power
    [~, inhSel] = max(peakPow(:,1:2), [], 2);  % 1 or 2
    inhaleHz = peakHz(sub2ind(size(peakHz), (1:size(peakHz,1))', inhSel));

    % exhale: choose the frequency from epochs 3:4 with larger peak power
    [~, exSelLocal] = max(peakPow(:,3:4), [], 2); % 1 or 2
    exSel = exSelLocal + 2;                       % 3 or 4
    exhaleHz = peakHz(sub2ind(size(peakHz), (1:size(peakHz,1))', exSel));

    inhaleExhaleHz = [inhaleHz exhaleHz];

    % subject-level 0-centering using mean across all final inhale/exhale points
    subjectGrandMeanHz = mean(inhaleExhaleHz, 'all', 'omitnan');
    inhaleExhaleHz = inhaleExhaleHz - subjectGrandMeanHz;

    % summary info (on centered values)
    subInfo.nBreaths(i)          = size(inhaleExhaleHz,1);
    subInfo.subjectGrandMeanHz(i)= subjectGrandMeanHz;
    subInfo.meanInhaleHz(i)      = mean(inhaleExhaleHz(:,1), 'omitnan');
    subInfo.sdInhaleHz(i)        = std(inhaleExhaleHz(:,1), 'omitnan');
    subInfo.meanExhaleHz(i)      = mean(inhaleExhaleHz(:,2), 'omitnan');
    subInfo.sdExhaleHz(i)        = std(inhaleExhaleHz(:,2), 'omitnan');

    thisCol = dotColors(i,:);
    lineCol = thisCol;

    % one jitter offset per breath, shared across both plotted points
    jit = (rand(size(inhaleExhaleHz,1),1) - 0.5) * 0.22;

    % breath-wise lines and dots
    for r = 1:size(inhaleExhaleHz,1)
        xPlot = xBase + jit(r);
        yPlot = inhaleExhaleHz(r,:);

        plot(ax, xPlot, yPlot, '-', 'Color', [lineCol, .1], 'LineWidth', 3, ...
            'HandleVisibility', 'off');

        scatter(ax, xPlot, yPlot, 28, ...
            'MarkerFaceColor', thisCol, ...
            'MarkerEdgeColor', thisCol, ...
            'MarkerFaceAlpha', 0.1, ...
            'MarkerEdgeAlpha', 0.1, ...
            'HandleVisibility', 'off');
    end

    % thicker subject mean line
    meanPlot = [subInfo.meanInhaleHz(i), subInfo.meanExhaleHz(i)];
    plot(ax, xBase, meanPlot, '-', ...
        'Color', thisCol, ...
        'LineWidth', 8, ...
        'HandleVisibility', 'off');

    % mean markers
    scatter(ax, xBase, meanPlot, 90, ...
        'MarkerFaceColor', thisCol, ...
        'MarkerEdgeColor', labCol, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    allPeakVals = [allPeakVals; inhaleExhaleHz(:); meanPlot(:)];
end

% drop empty summary rows if any cells were skipped/empty
subInfo = subInfo(subInfo.nBreaths > 0, :);

% --- axis formatting ---
xlim([0.5 2.5])
xticks(1:2)
xticklabels(phaseNames)

ylabel('Centered peak frequency (Hz)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

xlabel('Respiratory epoch', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

% title('Breath-wise peak frequency by respiratory epoch (subject-centered)', ...
%     'FontSize', 20, ...
%     'FontWeight', 'bold', ...
%     'FontName', 'Dotum', ...
%     'Color', labCol);

set(ax, 'Layer', 'top')
ylim([-7 10])
% sensible y-limits
% if ~isempty(allPeakVals)
%     ylim([floor(min(allPeakVals)-1), ceil(max(allPeakVals)+1)])
% else
%     ylim([-5 5])
% end

end