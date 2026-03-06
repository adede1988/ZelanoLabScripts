function hFig = plotFlatSpec_bySubjectBlocks(idxVec, subjectInfo, dataCell)
% Drop-in wrapper around the exact block-labeling logic.
% Keeps all hard-codes exactly as previously written.
frex = logspace(log10(.1),log10(200),300);
% --- pull cells + subject names for the selected rows ---
dup1Cells = dataCell(idxVec);
subNames  = string(subjectInfo(idxVec,8));
subNames  = strtrim(subNames);

% --- preprocess each cell exactly like your pipeline (mean over dim2, drop NaN rows) ---
cleanCells = cell(size(dup1Cells));
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

    cleanCells{i} = tmp;
end

% drop empties (in case any subject ends up with 0 trials after cleaning)
keep = ~cellfun(@isempty, cleanCells) & strlength(subNames) > 0;
cleanCells = cleanCells(keep);
subNames   = subNames(keep);

% --- group by subject name (stable order), concatenate within-subject ---
[uSub, ~, g] = unique(subNames, 'stable');

subMat = cell(numel(uSub),1);
subLen = zeros(numel(uSub),1);
for s = 1:numel(uSub)
    subMat{s} = cat(1, cleanCells{g==s});
    subLen(s) = size(subMat{s},1);
end

dup1FlatSpec = cat(1, subMat{:});  % [allTrials x nFreq]

% --- y-block bookkeeping ---
endsY   = cumsum(subLen);
startsY = [1; endsY(1:end-1)+1];
midsY   = startsY + (subLen-1)/2;

% --- plot ---
hFig = figure;
imagesc(dup1FlatSpec(:, frex>4));
plotFrex = frex(frex>4);

% (optional) make your xticklabels behave sensibly by setting xticks too
xt = 20:20:min(140, numel(plotFrex));
xticks(xt);
xticklabels(plotFrex(xt));

clim([-.5 1.5]);

% remove all y ticks/labels
set(gca, 'YTick', [], 'YTickLabel', []);

hold on;

% dashed separators between subjects
for s = 1:numel(uSub)-1
    yline(endsY(s)+0.5, 'k--', 'LineWidth', 1);
end

% subject labels centered in each block (placed just left of the image)
xText = 0.5;  % just left of first column (imagesc x starts at 1)
for s = 1:numel(uSub)
    text(xText, midsY(s), char(uSub(s)), ...
        'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
        'FontWeight','bold', 'Color','k', 'Interpreter','none', ...
        'Clipping','off');
end

hold off;
end