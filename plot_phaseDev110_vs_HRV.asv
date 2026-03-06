function hFig = plot_phaseDev110_vs_HRV(idxVec, allSubIDs, allBehDat)
% One-off helper for this pipeline.
% Scatter breath-wise phase deviation from 110° vs HRV index, colored by subject/session.
% Also overlays a best-fit line per subject/session.

% ---- hard codes ----
refDeg = 110;

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

for s = 1:nSub
    T = allBehDat{repRow(s)};

    % --- phase (deg), auto-detect radians ---
    phi = double(T.gamRspPACphase(:));
    mPhi = isfinite(phi);
    phi = phi(mPhi);

    if isempty(phi), continue; end

    if max(abs(phi)) <= (2*pi + 0.25)
        phiDeg = phi * 180/pi;
    else
        phiDeg = phi;
    end

    % signed circular deviation from 110° in [-180, 180]
    dDeg = atan2d(sind(phiDeg - refDeg), cosd(phiDeg - refDeg));

    % --- HRV index ---
    hrv = 1000*double(T.HRV_RMSSD30(:)) + 100*double(T.HRV_SDNN30(:)) + double(T.HRV_RSAamp(:));
    hrv = hrv(mPhi);  % match the finite phase mask applied above

    % remove NaNs/Infs jointly
    m = isfinite(dDeg) & isfinite(hrv);
    dDeg = dDeg(m);
    hrv  = hrv(m);

    if isempty(dDeg), continue; end

    % scatter
    scatter(ax, dDeg, hrv, 14, ...
        'MarkerFaceColor', cols(s,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.45);

    % % best-fit line (per subject)
    % if numel(dDeg) >= 2
    %     p = polyfit(dDeg, hrv, 1);
    %     xFit = linspace(min(dDeg), max(dDeg), 50);
    %     yFit = polyval(p, xFit);
    %     plot(ax, xFit, yFit, '-', 'Color', cols(s,:), 'LineWidth', 2);
    % end
end

xlabel(ax, 'Phase deviation from 110° (deg; signed shortest distance)');
ylabel(ax, 'HRV index = 1000*HRV\_RMSSD30 + 100*HRV\_SDNN30 + HRV\_RSAamp');
title(ax, 'Breath-wise phase deviation vs HRV (colored by subject/session)');

% optional legend (can get big; comment out if you don’t want it)
legend(ax, cellstr(uName), 'Interpreter','none', 'Box','off', 'Location','eastoutside');

box(ax,'off');
grid(ax,'on');
end