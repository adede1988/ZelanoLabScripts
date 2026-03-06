function hFig = plot_task_vs_HRVindex(idxVec, allSubIDs, allBehDat)
% One-off helper for this pipeline.
% X: task (categorical labels)
% Y: HRV index = 1000*HRV_RMSSD30 + 100*HRV_SDNN30 + HRV_RSAamp
% One line per subject/session (allSubIDs(:,7)), dots per task, lines connect dots.

% --- rows to use ---
if islogical(idxVec)
    sel = find(idxVec);
else
    sel = idxVec(:);
end

% --- subject/session id (col 7), name (col 8), session number (col 5) ---
subIDs = nan(numel(sel),1);
for i = 1:numel(sel)
    v = allSubIDs{sel(i),7};
    if isnumeric(v); subIDs(i) = double(v);
    else;            subIDs(i) = str2double(string(v));
    end
end
subNames = strtrim(string(allSubIDs(sel,8)));

sessNum = nan(numel(sel),1);
for i = 1:numel(sel)
    v = allSubIDs{sel(i),5};
    if isnumeric(v); sessNum(i) = double(v);
    else;            sessNum(i) = str2double(string(v));
    end
end

% --- unique subject/sessions (rows are redundant; keep first) ---
[uID, ia] = unique(subIDs, 'stable');
repRow = sel(ia);
uName  = subNames(ia);
uSess  = sessNum(ia);
nSub   = numel(uID);

% --- gather all tasks + per-subject means ---
allTasks = string.empty(0,1);
subTaskNames = cell(nSub,1);
subTaskMeans = cell(nSub,1);

for s = 1:nSub
    T = allBehDat{repRow(s)};

    % filter breaths
    if ismember('useVec', T.Properties.VariableNames)
        m = (T.useVec == 1);
    else
        m = true(height(T),1);
    end

    % tasks
    task = T.task;
    if iscell(task);       task = string(task);
    elseif iscategorical(task); task = string(task);
    else;                 task = string(task);
    end
    task = strtrim(task(m));

    % HRV index
    hrv = 1000*double(T.HRV_RMSSD30(m)) + 100*double(T.HRV_SDNN30(m)) + double(T.HRV_RSAamp(m));

    % finite joint mask
    mf = (strlength(task) > 0) & isfinite(hrv);
    task = task(mf);
    hrv  = hrv(mf);

    if isempty(task)
        subTaskNames{s} = string.empty(0,1);
        subTaskMeans{s} = [];
        continue
    end

    uT = unique(task, 'stable');
    mu = nan(numel(uT),1);
    for k = 1:numel(uT)
        mu(k) = mean(hrv(task==uT(k)), 'omitnan');
    end

    subTaskNames{s} = uT;
    subTaskMeans{s} = mu;

    allTasks = unique([allTasks; uT], 'stable');
end

nTask = numel(allTasks);

% --- build matrix [nSub x nTask] ---
Y = nan(nSub, nTask);
for s = 1:nSub
    uT = subTaskNames{s};
    mu = subTaskMeans{s};
    for k = 1:numel(uT)
        j = find(allTasks == uT(k), 1, 'first');
        Y(s,j) = mu(k);
    end
end

% --- plot ---
hFig = figure('Color','w');
ax = axes('Parent', hFig); hold(ax,'on');

cols = lines(max(nSub,1));
xpos = 1:nTask;

hLeg = gobjects(nSub,1);
for s = 1:nSub
    % line + markers; NaNs break lines where task missing
    hLeg(s) = plot(ax, xpos, Y(s,:), '-o', ...
        'Color', cols(s,:), ...
        'LineWidth', 1.75, ...
        'MarkerFaceColor', cols(s,:), ...
        'MarkerEdgeColor', cols(s,:), ...
        'MarkerSize', 5);
end

set(ax, 'XTick', xpos, 'XTickLabel', cellstr(allTasks), 'XTickLabelRotation', 35);

ylabel(ax, 'HRV index = 1000*HRV\_RMSSD30 + 100*HRV\_SDNN30 + HRV\_RSAamp');
xlabel(ax, 'Task');
title(ax, 'Mean HRV index by task (lines connect within subject/session)', 'Interpreter','none');

box(ax,'off');
grid(ax,'on');

% legend: subject name + session
legStr = strings(nSub,1);
for s = 1:nSub
    legStr(s) = sprintf('%s (sess%d) | ID %d', char(uName(s)), uSess(s), uID(s));
end
legend(ax, hLeg, cellstr(legStr), 'Interpreter','none', 'Box','off', 'Location','eastoutside');

end