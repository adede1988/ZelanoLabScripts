function hFig = plot_packStrength_vs_HRV(idxVec, allSubIDs, allBehDat)
% One-off helper for this pipeline.
% Scatter breath-wise gamRspPACKStrength vs HRV index, colored by subject/session.
% Overlays a best-fit line per subject/session.
%
% HRV index = 1000*HRV_RMSSD30 + 100*HRV_SDNN30 + HRV_RSAamp

% ---- resolve selected rows ----
if islogical(idxVec)
    sel = find(idxVec);
else
    sel = idxVec(:);
end

% ---- group by subject/session ID (allSubIDs col 7) ----
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

% representative row per subject/session (assumes behDat same across channels)
repRow = sel(ia);

cols = lines(nSub);

% ---- figure ----
hFig = figure('Color','w');
ax = axes('Parent', hFig); hold(ax,'on');

% store handles for a clean legend (one per subject)
hLeg = gobjects(nSub,1);

for s = 1:nSub
    T = allBehDat{repRow(s)};

    % --- PAC strength (raw) ---
    pac = double(T.gamRspPACLock(:));   %%%%%%%%%%%%%%%edit here%%%%%%%%%%%%%%%%%%%%%%%%%%

    % --- HRV index ---
    hrv = 1000*double(T.HRV_RMSSD30(:)) + 100*double(T.HRV_SDNN30(:)) + double(T.HRV_RSAamp(:));

    % joint finite mask
    m = isfinite(pac) & isfinite(hrv);
    pac = pac(m);
    hrv = hrv(m);

    if isempty(pac), continue; end

    % scatter
    hLeg(s) = scatter(ax, pac, hrv, 14, ...
        'MarkerFaceColor', cols(s,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.45);

    % best-fit line per subject
    if numel(pac) >= 2
        p = polyfit(pac, hrv, 1);
        xFit = linspace(min(pac), max(pac), 50);
        yFit = polyval(p, xFit);
        plot(ax, xFit, yFit, '-', 'Color', cols(s,:), 'LineWidth', 2, ...
            'HandleVisibility','off');
    end
end

xlabel(ax, 'gamRspPACKStrength (raw)');
ylabel(ax, 'HRV index = 1000*HRV\_RMSSD30 + 100*HRV\_SDNN30 + HRV\_RSAamp');
title(ax, 'Breath-wise PAC strength vs HRV (colored by subject/session)');
box(ax,'off');
grid(ax,'on');

% legend (skip empty handles if some subjects had no data)
mLeg = isgraphics(hLeg);
if any(mLeg)
    legend(ax, hLeg(mLeg), cellstr(uName(mLeg)), ...
        'Interpreter','none', 'Box','off', 'Location','eastoutside');
end
end