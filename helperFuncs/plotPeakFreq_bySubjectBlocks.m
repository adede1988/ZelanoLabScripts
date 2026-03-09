function [hFig, subInfo] = plotPeakFreq_bySubjectBlocks(idxVec, subjectInfo, dataCell, ...
    bgCol, labCol, dotColors)
% plotPeakFreq_bySubjectBlocks
%
% For each selected cell in dataCell:
%   - assumes data are [breath x block(5) x frequency]
%   - finds row-wise peak frequency in 25-60 Hz for each of the 5 block levels
%   - plots jittered dots with faint within-breath linking lines
%
% Colors indicate which dup1Cells index the breath came from.
%
% Optional:
%   dotColors : [nSelected x 3] RGB colors, one per selected cell

frex = logspace(log10(.1),log10(200),300);

phaseNames = {'inhale rise', 'inhale fall', 'exhale rise', 'exhale fall', 'pause'};
xBase = 1:5;

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
    'VariableNames', {'cellIdx','subName','nBreaths','meanPeakHz','sdPeakHz'});

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

    % row-wise peak frequency for each block level
    for b = 1:5
        cur = squeeze(tmp(:,b,:));   % [breath x nBand]
        if isvector(cur)
            cur = cur(:)';
        end

        goodRow = any(isfinite(cur), 2);
        if any(goodRow)
            [~, pkIdx] = max(cur(goodRow,:), [], 2);
            peakHz(goodRow,b) = fBand(pkIdx);
        end
    end

    % keep only breaths with all 5 block peaks defined
    keepRow = all(isfinite(peakHz), 2);
    peakHz  = peakHz(keepRow,:);

    if isempty(peakHz)
        continue
    end

    % summary info
    subInfo.nBreaths(i)   = size(peakHz,1);
    subInfo.meanPeakHz(i) = mean(peakHz, 'all', 'omitnan');
    subInfo.sdPeakHz(i)   = std(peakHz(:), 'omitnan');

    thisCol = dotColors(i,:);
    lineCol = 0.25*bgCol + 0.75*thisCol;   % faint line against dark bg
    lineCol = thisCol; 
    % one jitter offset per breath, shared across all 5 points for that breath
    jit = (rand(size(peakHz,1),1) - 0.5) * 0.22;

    for r = 1:size(peakHz,1)
        xPlot = xBase + jit(r);
        yPlot = peakHz(r,:);

        % faint linking line
        plot(ax, xPlot, yPlot, '-', 'Color', lineCol, 'LineWidth', 0.4, ...
            'HandleVisibility', 'off');

        % dots
        scatter(ax, xPlot, yPlot, 28, ...
            'MarkerFaceColor', thisCol, ...
            'MarkerEdgeColor', thisCol, ...
            'MarkerFaceAlpha', 0.6, ...
            'MarkerEdgeAlpha', 0.6, ...
            'HandleVisibility', 'off');
    end

    allPeakVals = [allPeakVals; peakHz(:)];
end

% drop empty summary rows if any cells were skipped/empty
subInfo = subInfo(subInfo.nBreaths > 0, :);

% --- axis formatting ---
xlim([0.5 5.5])
xticks(1:5)
xticklabels(phaseNames)

ylabel('Peak frequency (Hz)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

xlabel('Respiratory epoch', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

title('Breath-wise peak frequency by respiratory epoch', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

set(ax, 'Layer', 'top')

% sensible y-limits
if ~isempty(allPeakVals)
    ylim([max(24, floor(min(allPeakVals)-1)), min(61, ceil(max(allPeakVals)+1))])
else
    ylim([25 60])
end

end