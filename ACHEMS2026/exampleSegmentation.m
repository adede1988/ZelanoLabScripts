
% Run from debug of breathPiecewiseTemplateIdx.m
%assumes run through line 178
figSaveDir = 'G:\My Drive\cZelano\ACHEMS_2026\figs';
% --- styling/colors to match the other figures ---
bgCol    = [255 255 255]/255;
labCol   = [78 42 132]/255;
rspCol   = [195 176 163]/255;
xlineCol = [203 157 6]/255;

% --- derive final return index if it is not already in workspace ---
if ~exist('returnIdx','var')
    returnIdx = troughIdx + ca - 1;
end

% --- choose plotting window ---
wPlot = onsetIdx:(winEnd-1);
xPlot = x(wPlot);

% Handle either case:
%   1) t is full-length and matches x
%   2) t is already the windowed time vector for x(onsetIdx:winEnd-1)
if numel(t) == numel(x)
    tPlot = t(wPlot);

    t1 = t(peakIdx);
    t2 = t(crossBelowIdx);
    t3 = t(troughIdx);
    t4 = t(returnIdx);

elseif numel(t) == numel(wPlot)
    tPlot = t(:)';

    peakLocal       = peakIdx       - onsetIdx + 1;
    crossBelowLocal = crossBelowIdx - onsetIdx + 1;
    troughLocal     = troughIdx     - onsetIdx + 1;
    returnLocal     = returnIdx     - onsetIdx + 1;

    t1 = tPlot(peakLocal);
    t2 = tPlot(crossBelowLocal);
    t3 = tPlot(troughLocal);
    t4 = tPlot(returnLocal);

else
    error('t does not match either x or the plotted breath window.');
end

% --- figure/axes ---
hFig = figure( ...
    'Color', bgCol, ...
    'Position', [100 100 1100 500], ...
    'InvertHardcopy', 'off');

ax = axes('Parent', hFig);
hold(ax,'on')

ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.YColor     = labCol;
ax.TickDir    = 'out';
ax.TickLength = [0.018 0.018];
box(ax,'off')

% --- respiration trace ---
plot(ax, tPlot, xPlot, ...
    'Color', rspCol, ...
    'LineWidth', 8);

% --- epoch boundary lines ---
xline(ax, t1, '-', 'Color', xlineCol, 'LineWidth', 5, 'HandleVisibility', 'off');
xline(ax, t2, '-', 'Color', xlineCol, 'LineWidth', 5, 'HandleVisibility', 'off');
xline(ax, t3, '-', 'Color', xlineCol, 'LineWidth', 5, 'HandleVisibility', 'off');
xline(ax, t4, '-', 'Color', xlineCol, 'LineWidth', 5, 'HandleVisibility', 'off');

% --- labels/titles ---
xlabel(ax, 'Time (s)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

ylabel(ax, 'Respiration (flow rate au)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

% title(ax, 'Example respiratory cycle segmentation', ...
%     'FontSize', 20, ...
%     'FontWeight', 'bold', ...
%     'FontName', 'Dotum', ...
%     'Color', labCol);

% --- nice limits ---
xlim(ax, [tPlot(1) tPlot(end)])

yl = ylim(ax);
yr = range(yl);
ylim(ax, [yl(1) - 0.02*yr, yl(2) + 0.05*yr])

yl = ylim(ax);
yr = range(yl);

% place epoch labels just above the x-axis, inside the plotting area
textY = yl(1) + 0.05*yr;

xStart = tPlot(1);
xEnd   = tPlot(end);

epochMids = [ ...
    (xStart + t1)/2
    (t1 + t2)/2
    (t2 + t3)/2
    (t3 + t4)/2
    (t4 + xEnd)/2 ];

epochLabs = { ...
    'inhale rise'
    'inhale fall'
    'exhale rise'
    'exhale fall'
    'pause' };

for ii = 1:numel(epochMids)
    text(ax, epochMids(ii), textY, epochLabs{ii}, ...
        'Color', labCol, ...
        'FontSize', 13, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Clipping', 'on');
end

set(ax, 'Layer', 'top')

outFile = fullfile(figSaveDir, 'segmentationExample.png');


exportgraphics(hFig, outFile,'BackgroundColor', bgCol, 'Resolution', 300);
