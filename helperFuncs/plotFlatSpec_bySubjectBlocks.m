function [hFig, subInfo] = plotFlatSpec_bySubjectBlocks(idxVec, subjectInfo, dataCell, ...
    bgCol, labCol, labelColors)
% Optional input:
%   labelColors = [N x 3] RGB matrix, one row per displayed subject block.
%   If missing or wrong size, labCol is used for all subject labels.

% Drop-in wrapper around the exact block-labeling logic.
frex = logspace(log10(.1),log10(200),300);

% --- pull cells + subject names for the selected rows ---
dup1Cells = dataCell(idxVec);
subNames  = string(subjectInfo(idxVec,8));
subNames  = strtrim(subNames);
subCodes  = string(subjectInfo(idxVec,9));
subCodes  = strtrim(subCodes);

% --- preprocess each cell exactly like your pipeline (mean over dim2, drop NaN rows) ---
cleanCells = cell(size(dup1Cells));

% preallocate subject-level summary table
subInfo = table( ...
    strings(numel(dup1Cells),1), ...
    nan(numel(dup1Cells),1), ...
    nan(numel(dup1Cells),1), ...
    'VariableNames', {'subName','peakProm','peakSD'});

for i = 1:numel(dup1Cells)
    tmp = dup1Cells{i};

    % match your pipeline: average across dim 2, then squeeze
    tmp = squeeze(mean(tmp, 2, 'omitnan'));

    % match your pipeline: drop trials with NaN in column 3
    % (fallback: if <3 cols, drop any-NaN rows)
    if size(tmp,2) >= 3
        tmp(isnan(tmp(:,3)), :) = [];
    else
        tmp(any(isnan(tmp),2), :) = [];
    end

    test = max(tmp(:,219:252), [], 2);

    curName = string(subNames(i));
    curMean = mean(test, 'omitnan');
    curSD   = std(test,  'omitnan');

    disp(['subName: ' char(curName) ...
          ' peakProm: ' num2str(curMean) ...
          ' peakSD: ' num2str(curSD)])

    subInfo.subName(i)  = curName;
    subInfo.peakProm(i) = curMean;
    subInfo.peakSD(i)   = curSD;

    cleanCells{i} = tmp;
end

% drop empties (in case any subject ends up with 0 trials after cleaning)
keep  = ~cellfun(@isempty, cleanCells);
keep2 = strlength(subNames) > 0;
keep  = keep(:) & keep2(:);

cleanCells = cleanCells(keep);
subNames   = subNames(keep);
subCodes   = subCodes(keep);

% keep subInfo aligned with retained rows
subInfo = subInfo(keep,:);

% --- group by subject name (stable order), concatenate within-subject ---
[uSub, ~, g] = unique(subNames, 'stable');

subMat = cell(numel(uSub),1);
subLen = zeros(numel(uSub),1);
uCod   = strings(numel(uSub),1);

for s = 1:numel(uSub)
    subMat{s} = cat(1, cleanCells{g==s});
    subLen(s) = size(subMat{s},1);

    % grab the first code associated with this subject block
    firstIdx = find(g==s, 1, 'first');
    uCod(s) = subCodes(firstIdx);
end

dup1FlatSpec = cat(1, subMat{:});  % [allTrials x nFreq]

% --- y-block bookkeeping ---
endsY   = cumsum(subLen);
startsY = [1; endsY(1:end-1)+1];
midsY   = startsY + (subLen-1)/2;

% --- label color handling ---
if nargin < 6 || isempty(labelColors) || ...
        ~isnumeric(labelColors) || size(labelColors,2) ~= 3 || ...
        size(labelColors,1) ~= numel(uSub)
    useLabelColors = repmat(labCol, numel(uSub), 1);
else
    useLabelColors = labelColors;
end

% --- plot ---
hFig = figure('position', [0,0,600, 370]);
ax = axes;
ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.YColor     = labCol;
box(ax, 'off');

imagesc(dup1FlatSpec(:, frex>4));
plotFrex = frex(frex>4);

xt = 20:20:min(140, numel(plotFrex));
xticks(xt)
xticklabels(string(round(plotFrex(xt))))
set(gca, 'FontSize', 20, 'FontWeight', 'bold', ...
    'FontName', 'Dotum', 'XColor', labCol)

clim([-.5 1.5]);

% remove all y ticks/labels
set(gca, 'YTick', [], 'YTickLabel', []);

hold on;

% dashed separators between subjects
for s = 1:numel(uSub)-1
    yline(endsY(s)+0.5, 'k--', 'LineWidth', 1);
end

% subject labels centered in each block
xText = 0.5;
for s = 1:numel(uSub)
    text(xText, midsY(s), char(uCod(s)), ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'Interpreter','none', ...
        'Clipping','off', ...
        'FontSize', 20, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', useLabelColors(s,:));
end