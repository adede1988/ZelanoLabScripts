function [hFig, subInfo, diffCell, ampCell] = ...
    plotInhaleExhaleFreqDiff_vs_BreathAmp(idxVec, subjectInfo, dataCell, ...
    allBehDat, allResp, bgCol, labCol, dotColors)
% plotInhaleExhaleFreqDiff_vs_BreathAmp
%
% For each selected cell in dataCell:
%   - assumes spectral data are [breath x block(5) x frequency]
%   - computes, for each breath, the inhale-exhale frequency difference:
%         inhale = frequency from whichever of epochs 1:2 has larger peak power
%         exhale = frequency from whichever of epochs 3:4 has larger peak power
%         diffHz = inhaleHz - exhaleHz
%   - computes breath amplitude from respiratory flow:
%         amplitude = max(flow during breath) - min(flow during breath)
%     where:
%         onset sample = 1000
%         end sample   = 1000 + length*500
%   - makes a scatter plot of:
%         x = inhale-exhale frequency difference
%         y = breath amplitude
%     with dots colored by participant
%
% Inputs:
%   idxVec      : indices into subjectInfo / dataCell / allBehDat / allResp
%   subjectInfo : metadata array/table; subject name assumed in column 8
%   dataCell    : cell array, each cell = [breath x 5 x frequency]
%   allBehDat   : cell array of tables, one per entry; must contain column 'length'
%   allResp     : cell array of respiration matrices [breath x 6000]
%   bgCol       : background color
%   labCol      : axis/label color
%   dotColors   : optional [nSelected x 3] RGB colors
%
% Outputs:
%   hFig     : figure handle
%   subInfo  : summary table, one row per retained selected entry
%   diffCell : per-entry cell array of inhale-exhale frequency differences
%   ampCell  : per-entry cell array of breath amplitudes
%
% Notes:
%   - Sampling rate is assumed to be 500 Hz
%   - Breath onset is assumed to be at sample 1000
%   - Breath end is computed as round(1000 + length*500)
%   - If multiple selected rows belong to the same subject name, they are
%     plotted with the same color (taken from the first occurrence)

frex = logspace(log10(.1), log10(200), 300);
fs = 500;
breathOnsetSamp = 1000;
nRespSamp = 6000;

% --- pull selected cells + names ---
dup1Cells = dataCell(idxVec);
behCells  = allBehDat(idxVec);
respCells = allResp(idxVec);

subNames = string(subjectInfo(idxVec,8));
subNames = strtrim(subNames);

nCell = numel(dup1Cells);

% optional colors
if nargin < 8 || isempty(dotColors)
    dotColors = lines(nCell);
elseif size(dotColors,1) < nCell
    error('dotColors must have at least as many rows as selected cells.');
end

% force repeated subject names to use same color
uNames = unique(subNames, 'stable');
for u = 1:numel(uNames)
    rows = find(subNames == uNames(u));
    if numel(rows) > 1
        dotColors(rows,:) = repmat(dotColors(rows(1),:), numel(rows), 1);
    end
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
    'VariableNames', {'cellIdx','subName','nBreaths','meanDiffHz','sdDiffHz','meanAmp','sdAmp'});

diffCell = cell(nCell,1);
ampCell  = cell(nCell,1);

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

allX = [];
allY = [];

% --- loop over selected entries ---
for i = 1:nCell
    tmp = dup1Cells{i};
    beh = behCells{i};
    rsp = respCells{i};

    % basic checks
    if ndims(tmp) ~= 3
        warning('Skipping cell %d: expected spectral array [breath x block x freq].', i);
        continue
    end
    if size(tmp,2) < 5
        warning('Skipping cell %d: second dimension has fewer than 5 levels.', i);
        continue
    end
    if ~istable(beh) || ~ismember('length', beh.Properties.VariableNames)
        warning('Skipping cell %d: allBehDat entry must be a table with a column named length.', i);
        continue
    end
    if ~ismatrix(rsp) || size(rsp,2) ~= nRespSamp
        warning('Skipping cell %d: allResp entry must be [breath x 6000].', i);
        continue
    end

    nBreathSpec = size(tmp,1);
    nBreathBeh  = height(beh);
    nBreathResp = size(rsp,1);

    if nBreathSpec ~= nBreathBeh || nBreathSpec ~= nBreathResp
        warning('Skipping cell %d: number of breaths does not match across dataCell/allBehDat/allResp.', i);
        continue
    end

    % keep only first 5 levels of dim 2, and 25-60 Hz in dim 3
    tmp = tmp(:,1:5,fMask);   % [breath x 5 x nBand]

    peakHz  = nan(nBreathSpec, 5);
    peakPow = nan(nBreathSpec, 5);

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

    % inhale: choose frequency from epochs 1:2 with larger peak power
    inhaleHz = nan(nBreathSpec,1);
    exhaleHz = nan(nBreathSpec,1);

    good12 = all(isfinite(peakPow(:,1:2)), 2) & all(isfinite(peakHz(:,1:2)), 2);
    if any(good12)
        [~, inhSel] = max(peakPow(good12,1:2), [], 2);   % 1 or 2
        tmpHz = peakHz(good12,1:2);
        inhaleHz(good12) = tmpHz(sub2ind(size(tmpHz), (1:size(tmpHz,1))', inhSel));
    end

    good34 = all(isfinite(peakPow(:,3:4)), 2) & all(isfinite(peakHz(:,3:4)), 2);
    if any(good34)
        [~, exSel] = max(peakPow(good34,3:4), [], 2);    % 1 or 2 within local block
        tmpHz = peakHz(good34,3:4);
        exhaleHz(good34) = tmpHz(sub2ind(size(tmpHz), (1:size(tmpHz,1))', exSel));
    end

    diffHz = inhaleHz - exhaleHz;

    % respiratory amplitude during the breath window
    breathAmp = nan(nBreathSpec,1);
    breathLenSec = beh.length;

    for r = 1:nBreathSpec
        L = breathLenSec(r);

        if ~isfinite(L) || L <= 0
            continue
        end

        endSamp = round(breathOnsetSamp + L*fs);
        endSamp = min(max(endSamp, breathOnsetSamp), nRespSamp);

        seg = rsp(r, breathOnsetSamp:endSamp);
        seg = seg(isfinite(seg));

        if isempty(seg)
            continue
        end

        breathAmp(r) = max(seg) - min(seg);
    end

    keepRow = isfinite(diffHz) & isfinite(breathAmp);

    diffHz = diffHz(keepRow);
    breathAmp = breathAmp(keepRow);

    if isempty(diffHz)
        continue
    end

    diffCell{i} = diffHz;
    ampCell{i}  = breathAmp;

    subInfo.nBreaths(i)   = numel(diffHz);
    subInfo.meanDiffHz(i) = mean(diffHz, 'omitnan');
    subInfo.sdDiffHz(i)   = std(diffHz, 'omitnan');
    subInfo.meanAmp(i)    = mean(breathAmp, 'omitnan');
    subInfo.sdAmp(i)      = std(breathAmp, 'omitnan');

    thisCol = dotColors(i,:);

    scatter(ax, diffHz, breathAmp, 150, ...
        'MarkerFaceColor', thisCol, ...
        'MarkerEdgeColor', thisCol, ...
        'MarkerFaceAlpha', 0.55, ...
        'MarkerEdgeAlpha', 0.55, ...
        'HandleVisibility', 'off');

    % thicker subject mean point
    scatter(ax, subInfo.meanDiffHz(i), subInfo.meanAmp(i), 600, ...
        'MarkerFaceColor', thisCol, ...
        'MarkerEdgeColor', labCol, ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');

    allX = [allX; diffHz(:); subInfo.meanDiffHz(i)];
    allY = [allY; breathAmp(:); subInfo.meanAmp(i)];
end

% drop empty rows
keepSub = subInfo.nBreaths > 0;
subInfo = subInfo(keepSub,:);
diffCell = diffCell(keepSub);
ampCell  = ampCell(keepSub);

xlabel('Inhale - exhale peak frequency difference (Hz)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

ylabel('Breath amplitude (max - min flow)', ...
    'FontSize', 20, ...
    'FontWeight', 'bold', ...
    'FontName', 'Dotum', ...
    'Color', labCol);

% title('Breath amplitude vs inhale - exhale frequency difference', ...
%     'FontSize', 20, ...
%     'FontWeight', 'bold', ...
%     'FontName', 'Dotum', ...
%     'Color', labCol);

set(ax, 'Layer', 'top')

% sensible limits
if ~isempty(allX)
    xPad = max(0.5, 0.05 * range(allX));
    if ~isfinite(xPad), xPad = 1; end
    xlim([min(allX)-xPad, max(allX)+xPad]);
end

if ~isempty(allY)
    yPad = max(0.05 * range(allY), eps);
    if ~isfinite(yPad), yPad = 1; end
    ylim([max(0, min(allY)-yPad), max(allY)+yPad]);
end

end