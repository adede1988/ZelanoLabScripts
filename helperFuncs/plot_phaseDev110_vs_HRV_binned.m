function hFig = plot_phaseDev110_vs_HRV_binned(idxVec, allSubIDs, allBehDat)
% One-off helper for this pipeline.
% For each subject/session: bin phase deviation from 110° (signed, circular) into
% symmetric bins around 0 (center = ±5°, then ±(5..15), ±(15..25), ...),
% compute mean HRV index per bin, and plot a smooth(ish) line through bin-centers.
%
% One line per subject/session (colored). No scatter.

% ---- hard codes ----
refDeg  = 110;
binHalf = 5;          % center bin is [-5,+5]
binW    = 30;         % flanking bins are 10° wide: [5..15], [15..25], ...
maxAbs  = 175;        % stop before 180 to avoid edge weirdness

% bin centers and edges
centers = (-maxAbs:binW:maxAbs).';
edges   = [centers - binW/2, centers + binW/2]; % [nBin x 2]
% force center bin to be exactly [-5,+5]
[~, i0] = min(abs(centers));
edges(i0,:) = [-binHalf, binHalf];
centers(i0) = 0;

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

[uID, ia] = unique(subIDs, 'stable');
uName = subNames(ia);
nSub  = numel(uID);

% representative row per subject/session (assumes behDat same across channels)
repRow = sel(ia);

cols = lines(nSub);

% ---- figure ----
hFig = figure('Color','w');
ax = axes('Parent', hFig); hold(ax,'on');

hLeg = gobjects(nSub,1);

for s = 1:nSub
    T = allBehDat{repRow(s)};

    % useVec==1 only (per your pipeline convention)
    if ismember('useVec', T.Properties.VariableNames)
        mUse = (T.useVec == 1);
    else
        mUse = true(height(T),1);
    end

    % --- phase (deg), auto-detect radians ---
    phi = double(T.gamRspPACphase(:));
    mPhi = isfinite(phi) & mUse;
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
    hrv = hrv(mPhi);  % match mask

    % joint finite mask
    m = isfinite(dDeg) & isfinite(hrv);
    dDeg = dDeg(m);
    hrv  = hrv(m);

    if isempty(dDeg), continue; end

    % --- bin means ---
    yBin = nan(numel(centers),1);
    for b = 1:numel(centers)
        lo = edges(b,1);
        hi = edges(b,2);

        if b == i0
            mb = (dDeg >= lo) & (dDeg <= hi);       % include edges for center
        else
            mb = (dDeg >  lo) & (dDeg <= hi);       % half-open elsewhere
        end

        if any(mb)
            yBin(b) = mean(hrv(mb), 'omitnan');
        end
    end

    % remove empty bins
    mB = isfinite(yBin);
    if nnz(mB) < 2, continue; end

    xPlot = centers(mB);
    yPlot = yBin(mB);

    % smooth-ish line: interpolate onto a dense grid then movmean
    xFine = linspace(min(xPlot), max(xPlot), 200);
    yFine = interp1(xPlot, yPlot, xFine, 'linear', 'extrap');

    % moving average smoothing across the dense grid
    win = 11; % odd window
    yFine = movmean(yFine, win);

    hLeg(s) = plot(ax, xFine, yFine, '-', 'Color', cols(s,:), 'LineWidth', 2);
end

xlabel(ax, 'Phase deviation from 110° (deg; signed shortest distance)');
ylabel(ax, 'HRV index = 1000*HRV\_RMSSD30 + 100*HRV\_SDNN30 + HRV\_RSAamp');
title(ax, 'Binned mean HRV vs phase deviation (one line per subject/session)', 'Interpreter','none');

% legend (skip empty handles)
mH = isgraphics(hLeg);
if any(mH)
    legend(ax, hLeg(mH), cellstr(uName(mH)), 'Interpreter','none', 'Box','off', 'Location','eastoutside');
end

box(ax,'off');
grid(ax,'on');
xline(ax, 0, 'k:', 'LineWidth', 1);  % target alignment reference
end