function [hFig, subTab] = plot_gamRspPACphase_circDensity(idxVec, allSubIDs, allBehDat, boldSubID)
% Drop-in helper (NOT general purpose).
% Plots lightly-shaded gray circular density regions for each subject (gamRspPACphase),
% with colored resultant vectors on top. Optionally "bold" one subject's density region
% by passing their numeric subject/session ID (allSubIDs(:,7)).
%
% Also returns a subjects x 10 table:
%   - all 8 columns of allSubIDs for the representative row of each subject/session
%   - resultant vector length (R)
%   - resultant vector direction (muDeg; 0..360 degrees)

if nargin < 4
    boldSubID = [];
end

% --- resolve selected rows ---
if islogical(idxVec)
    sel = find(idxVec);
else
    sel = idxVec(:);
end

% --- subject/session IDs (col 7) + names (col 8) ---
subIDs = nan(numel(sel),1);
for i = 1:numel(sel)
    v = allSubIDs{sel(i),7};
    if isnumeric(v); subIDs(i) = double(v);
    else;            subIDs(i) = str2double(string(v));
    end
end
subNames = strtrim(string(allSubIDs(sel,8)));

[uID, ia, g] = unique(subIDs, 'stable');   % unique subject/session IDs
uName = subNames(ia);
nSub  = numel(uID);

% representative row in allSubIDs / allBehDat for each subject/session
repRow = sel(ia);  % absolute row index

% --- circular density settings (von Mises KDE) ---
theta = linspace(-pi, pi, 361);
kappa = 8;
normC = 2*pi*besseli(0,kappa);

% colors for vectors
cols = lines(nSub);

% which subject to bold
boldIdx = [];
if ~isempty(boldSubID)
    boldIdx = find(uID == double(boldSubID), 1, 'first');
end

% --- precompute phase vectors per subject/session (so plot + table match) ---
phiCell = cell(nSub,1);
for s = 1:nSub
    T = allBehDat{repRow(s)};
    phi = double(T.gamRspPACphase(:));
    phi = phi(isfinite(phi));
    if isempty(phi)
        phiCell{s} = [];
        continue
    end

    % auto-detect degrees vs radians
    if max(abs(phi)) > (2*pi + 0.25)
        phi = phi * pi/180;
    end

    % wrap to [-pi, pi]
    phiCell{s} = atan2(sin(phi), cos(phi));
end

% --- figure / axes (Cartesian; we draw polar geometry ourselves) ---
hFig = figure('Color','w');
ax = axes('Parent', hFig);
hold(ax,'on');
axis(ax,'equal');
axis(ax,'off');

% faint unit circle
tt = linspace(-pi, pi, 400);
plot(ax, cos(tt), sin(tt), 'k:', 'LineWidth', 1);

% ---- edge tick labels ----
tickDeg  = [0 90 180 270];
tickText = { ...
    '0 = inhale peak', ...
    '90 = inhale/exhale transition', ...
    '180 = exhale trough', ...
    '270 = exhale/inhale transition' ...
    };

rLab = 1.14;
rTck = 1.02;
rTck2= 0.98;

for k = 1:numel(tickDeg)
    ang = tickDeg(k) * pi/180;
    [x1,y1] = pol2cart(ang, rTck2);
    [x2,y2] = pol2cart(ang, rTck);
    plot(ax, [x1 x2], [y1 y2], 'k-', 'LineWidth', 1);

    [xl,yl] = pol2cart(ang, rLab);

    if tickDeg(k)==0
        ha = 'left';  va = 'middle';
    elseif tickDeg(k)==180
        ha = 'right'; va = 'middle';
    elseif tickDeg(k)==90
        ha = 'center'; va = 'bottom';
    else % 270
        ha = 'center'; va = 'top';
    end

    text(ax, xl, yl, tickText{k}, ...
        'HorizontalAlignment', ha, 'VerticalAlignment', va, ...
        'Interpreter','none', 'FontSize', 9, 'Color','k');
end

% --- draw all gray density patches first (so they accumulate) ---
grayCol   = [0.5 0.5 0.5];
grayAlpha = 0.05;

for s = 1:nSub
    phi = phiCell{s};
    if isempty(phi), continue; end

    dens = mean(exp(kappa*cos(theta(:) - phi(:)')), 2) ./ normC;
    dens = dens ./ max(dens);

    thetaPoly = [theta, fliplr(theta)];
    rPoly     = [dens(:).', zeros(1, numel(theta))];
    [xp, yp]  = pol2cart(thetaPoly, rPoly);

    p = patch(ax, xp, yp, grayCol, 'EdgeColor','none', 'FaceAlpha', grayAlpha);
    p.HandleVisibility = 'off';
end

% --- if requested, redraw the bold subject density patch in its color ---
if ~isempty(boldIdx)
    phi = phiCell{boldIdx};
    if ~isempty(phi)
        dens = mean(exp(kappa*cos(theta(:) - phi(:)')), 2) ./ normC;
        dens = dens ./ max(dens);

        thetaPoly = [theta, fliplr(theta)];
        rPoly     = [dens(:).', zeros(1, numel(theta))];
        [xp, yp]  = pol2cart(thetaPoly, rPoly);

        pb = patch(ax, xp, yp, cols(boldIdx,:), ...
            'EdgeColor', cols(boldIdx,:), 'LineWidth', 1.5, ...
            'FaceAlpha', 0.25);
        pb.HandleVisibility = 'off';
    end
end

% --- draw resultant vectors on top (colored, thicker) + store for table ---
hVec  = gobjects(nSub,1);     % for legend
hasVec= false(nSub,1);

Rvec    = nan(nSub,1);
muDeg   = nan(nSub,1);        % 0..360 degrees (matches tick semantics)

for s = 1:nSub
    phi = phiCell{s};
    if isempty(phi), continue; end

    m   = mean(exp(1i*phi));
    mu  = angle(m);           % radians
    R   = abs(m);

    Rvec(s)  = R;
    muDeg(s) = mod(mu * 180/pi, 360);

    [xv, yv] = pol2cart(mu, R);

    lw = 2.75;
    if ~isempty(boldIdx) && s == boldIdx
        lw = 3.5;
    end

    hVec(s) = plot(ax, [0 xv], [0 yv], '-', 'Color', cols(s,:), 'LineWidth', lw);
    plot(ax, xv, yv, 'o', 'MarkerSize', 5, ...
        'MarkerFaceColor', cols(s,:), 'MarkerEdgeColor', cols(s,:), ...
        'HandleVisibility','off');

    hasVec(s) = true;
end

% limits
xlim(ax, [-1.25 1.25]);
ylim(ax, [-1.20 1.20]);

% title
if ~isempty(boldIdx)
    title(ax, sprintf('gamRspPACphase densities (gray) + resultant vectors | bold: %s (ID %d)', ...
        char(uName(boldIdx)), uID(boldIdx)), 'Interpreter','none');
else
    title(ax, 'gamRspPACphase densities (gray) + resultant vectors (colored)', 'Interpreter','none');
end

% legend with subject names keyed to vector colors
if any(hasVec)
    legend(ax, hVec(hasVec), cellstr(uName(hasVec)), ...
        'Interpreter','none', 'Box','off', 'Location','eastoutside');
end

% --- subjects x 10 table output ---
repInfo = allSubIDs(repRow, 1:8);  % 8 columns per subject/session
subTab  = cell2table(repInfo, 'VariableNames', ...
    {'col1','col2','col3','col4','col5','col6','col7','col8'});

subTab.R     = Rvec;
subTab.muDeg = muDeg;

end