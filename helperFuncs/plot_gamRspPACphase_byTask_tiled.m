function hFig = plot_gamRspPACphase_byTask_tiled(idxVec, allSubIDs, allBehDat)
% One-off helper for this pipeline.
% Tiled layout: one panel per subject/session (unique allSubIDs(:,7) within idxVec).
% Within each panel: split breaths by T.task (task column inside the table) and plot:
%   - circular density (lightly shaded) per task
%   - resultant vector (thicker) per task
% Legend indicates task. Title indicates subject/session.

% ----- selection rows -----
if islogical(idxVec)
    sel = find(idxVec);
else
    sel = idxVec(:);
end

% ----- subject/session IDs (col 7) + name (col 8) + session number (col 5) -----
subIDs = nan(numel(sel),1);
for i = 1:numel(sel)
    v = allSubIDs{sel(i),7};
    if isnumeric(v); subIDs(i) = double(v);
    else;            subIDs(i) = str2double(string(v));
    end
end
subNames = strtrim(string(allSubIDs(sel,8)));

[uID, ia, g] = unique(subIDs, 'stable');
uName = subNames(ia);
nSub  = numel(uID);
repRow = sel(ia);  % one representative row per subject/session

% ----- density settings -----
theta = linspace(-pi, pi, 361);
kappa = 8;
normC = 2*pi*besseli(0,kappa);

% ----- tiled layout -----
nCol = ceil(sqrt(nSub));
nRow = ceil(nSub / nCol);

hFig = figure('Color','w');
tlo  = tiledlayout(hFig, nRow, nCol, 'TileSpacing','compact', 'Padding','compact');

% figure-level annotation for phase semantics
annotation(hFig,'textbox',[0 0.965 1 0.035], ...
    'String','Phase ticks: 0=inhale peak; 90=inhale/exhale transition; 180=exhale trough; 270=exhale/inhale transition', ...
    'EdgeColor','none', 'HorizontalAlignment','center', 'Interpreter','none', 'FontSize',9);

% ----- draw each subject/session panel -----
for s = 1:nSub
    ax = nexttile(tlo, s);
    hold(ax,'on');
    axis(ax,'equal');
    axis(ax,'off');

    % base unit circle
    tt = linspace(-pi, pi, 400);
    plot(ax, cos(tt), sin(tt), 'k:', 'LineWidth', 1);

    % simple tick marks + numeric labels
    tickDeg = [0 90 180 270];
    rTck1 = 0.98; rTck2 = 1.02; rLab = 1.12;
    for k = 1:numel(tickDeg)
        ang = tickDeg(k)*pi/180;
        [x1,y1] = pol2cart(ang, rTck1);
        [x2,y2] = pol2cart(ang, rTck2);
        plot(ax, [x1 x2], [y1 y2], 'k-', 'LineWidth', 1);

        [xl,yl] = pol2cart(ang, rLab);
        text(ax, xl, yl, sprintf('%d', tickDeg(k)), ...
            'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
            'Interpreter','none', 'FontSize',8, 'Color','k');
    end

    % table for this subject/session
    T = allBehDat{repRow(s)};

    % phase (gamRspPACphase)
    phi = double(T.gamRspPACphase(:));
    mPhi = isfinite(phi);

    % optional useVec filter if present
    if ismember('useVec', T.Properties.VariableNames)
        mUse = (T.useVec == 1);
    else
        mUse = true(size(phi));
    end
    mKeep = mPhi & mUse;

    phi = phi(mKeep);
    if isempty(phi)
        xlim(ax,[-1.25 1.25]); ylim(ax,[-1.20 1.20]);
        title(ax, sprintf('%s | ID %d', char(uName(s)), uID(s)), 'Interpreter','none');
        continue
    end

    % auto-detect degrees vs radians, convert to radians if needed
    if max(abs(phi)) > (2*pi + 0.25)
        phi = phi * pi/180;
    end
    phi = atan2(sin(phi), cos(phi)); % wrap [-pi, pi]

    % task column (assumed to exist as 'task'; fallback to 'taskID' if needed)
    if ismember('task', T.Properties.VariableNames)
        task = T.task;
    else
        task = T.taskID;
    end
    task = task(mKeep);

    % normalize task type
    if iscell(task)
        task = string(task);
    end
    if iscategorical(task)
        task = string(task);
    end
    if isstring(task)
        task = strtrim(task);
    end

    % unique tasks (stable)
    uTask = unique(task, 'stable');
    nTask = numel(uTask);
    cols  = lines(max(nTask,1));

    % legend handles (use the vector lines for legend)
    hLeg = gobjects(nTask,1);
    hasLeg = false(nTask,1);

    for k = 1:nTask
        mk = (task == uTask(k));
        phik = phi(mk);
        phik = phik(isfinite(phik));
        if isempty(phik), continue; end

        % KDE (von Mises) and normalize to max=1 for display
        dens = mean(exp(kappa*cos(theta(:) - phik(:)')), 2) ./ normC;
        dens = dens ./ max(dens);

        % filled region polygon
        thetaPoly = [theta, fliplr(theta)];
        rPoly     = [dens(:).', zeros(1, numel(theta))];
        [xp, yp]  = pol2cart(thetaPoly, rPoly);

        patch(ax, xp, yp, cols(k,:), ...
            'EdgeColor', cols(k,:), 'LineWidth', 1.0, ...
            'FaceAlpha', 0.12, 'HandleVisibility','off');

        % resultant vector (thicker)
        m  = mean(exp(1i*phik));
        mu = angle(m);
        R  = abs(m);
        [xv,yv] = pol2cart(mu, R);

        hLeg(k) = plot(ax, [0 xv], [0 yv], '-', 'Color', cols(k,:), 'LineWidth', 3.0);
        plot(ax, xv, yv, 'o', 'MarkerSize', 4, ...
            'MarkerFaceColor', cols(k,:), 'MarkerEdgeColor', cols(k,:), ...
            'HandleVisibility','off');

        hasLeg(k) = true;
    end

    % title: subject/session
    ss = allSubIDs{repRow(s),5};
    if isnumeric(ss); ssn = double(ss); else; ssn = str2double(string(ss)); end
    title(ax, sprintf('%s (sess%d) | ID %d', char(uName(s)), ssn, uID(s)), 'Interpreter','none');

    % limits
    xlim(ax, [-1.25 1.25]);
    ylim(ax, [-1.20 1.20]);

    % legend: task labels
    if any(hasLeg)
        legend(ax, hLeg(hasLeg), cellstr(uTask(hasLeg)), ...
            'Interpreter','none', 'Box','off', 'Location','best');
    end
end

end