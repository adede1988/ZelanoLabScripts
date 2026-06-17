function [hFig, subTab, countMat, binCtrDeg] = plot_gamRspPhase_circStackedHist( ...
    idxVec, allSubIDs, allBehDat, colors, labCol, bgCol, phaseSource, nBins)
% Breath-level stacked circular histogram by subject/session.
%
% New helper based on plot_gamRspPACphase_circDensity, but instead of
% subject-level circular density regions it plots a stacked circular
% histogram where each bar reflects the total number of breaths in that
% phase bin, and the stacked colors show each subject's contribution.
%
% INPUTS
%   idxVec       : logical or numeric index into allSubIDs / allBehDat
%   allSubIDs    : cell array metadata
%   allBehDat    : cell array of tables/structs, one per selected row
%   colors       : [nSub x 3] RGB colors (optional)
%   labCol       : label color
%   bgCol        : background color
%   phaseSource  : 'gamRspPACphase' (default) or 'gamPeakidx50'
%   nBins        : number of circular histogram bins (default = 24)
%
% OUTPUTS
%   hFig         : figure handle
%   subTab       : one row per subject/session with metadata + summary stats
%   countMat     : [nSub x nBins] breath counts per subject per angular bin
%   binCtrDeg    : bin centers in degrees (0..360)
%
% NOTES
%   - For phaseSource = 'gamRspPACphase', data are assumed to be radians
%     unless values exceed ~2*pi, in which case they are treated as degrees.
%   - For phaseSource = 'gamPeakidx50', values 1..50 are mapped onto
%     0..360 degrees via:
%           deg = (idx - 1) * (360/50)
%     so:
%           1  -> 0 deg   (inhale start)
%           11 -> 72 deg  (inhale peak)
%           21 -> 144 deg (inhale/exhale transition)
%           31 -> 216 deg (exhale peak)
%           41 -> 288 deg (exhale/pause transition)

if nargin < 7 || isempty(phaseSource)
    phaseSource = 'gamRspPACphase';
end
if nargin < 8 || isempty(nBins)
    nBins = 24;
end

% ---------- resolve selected rows ----------
if islogical(idxVec)
    sel = find(idxVec);
else
    sel = idxVec(:);
end

% ---------- subject/session IDs (col 7) + names (col 8) ----------
subIDs = nan(numel(sel),1);
for i = 1:numel(sel)
    v = allSubIDs{sel(i),7};
    if isnumeric(v)
        subIDs(i) = double(v);
    else
        subIDs(i) = str2double(string(v));
    end
end
subNames = strtrim(string(allSubIDs(sel,8)));

[uID, ia, ~] = unique(subIDs, 'stable');
uName = subNames(ia);
nSub  = numel(uID);

% representative row in allSubIDs / allBehDat for each subject/session
repRow = sel(ia);

% ---------- colors ----------
if nargin < 4 || isempty(colors)
    colors = [];
end

if size(colors,1) == nSub
    cols = colors;
else
    warmCols = [ ...
        255 140  66;
        255  92  92;
        255 196  61;
        255 160 122;
        223 144 109;
        255 214 102;
        255 120 102;
        255 170  95] / 255;
    cols = warmCols(mod(0:nSub-1, size(warmCols,1)) + 1, :);
end

% ---------- source-specific setup ----------
phaseSource = char(string(phaseSource));

switch lower(phaseSource)
    case lower('gamRspPACphase')
        tickDeg  = [0 90 180 270];
        tickText = { ...
            'inhale start', ...
            'in/ex trans', ...
            'exhale trough', ...
            'ex/pause trans'};
        srcLabel = 'gamRspPACphase';

    case lower('gamPeakidx50')
        tickIdx  = [1 11 21 31 41];
        tickDeg  = (tickIdx - 1) * (360/50);
        tickText = { ...
            'inhale start', ...
            'inhale peak', ...
            'in/ex trans', ...
            'exhale peak', ...
            'ex/pause trans'};
        srcLabel = 'gamPeakidx50';

    otherwise
        error('phaseSource must be ''gamRspPACphase'' or ''gamPeakidx50''.');
end

% ---------- extract phase vectors per subject/session ----------
phiCell   = cell(nSub,1);
nBreaths  = zeros(nSub,1);
Rvec      = nan(nSub,1);
muDeg     = nan(nSub,1);

for s = 1:nSub
    T = allBehDat{repRow(s)};

    switch lower(phaseSource)
        case lower('gamRspPACphase')
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

            % wrap to [0, 2*pi)
            phi = mod(phi, 2*pi);

        case lower('gamPeakidx50')
            phiIdx = double(T.gamPeakidx50(:));
            phi = double(T.gamBackground(:));
            
            phiIdx = phiIdx(isfinite(phi));
            phiIdx = phiIdx(phiIdx >= 1 & phiIdx <= 50);

            if isempty(phiIdx)
                phiCell{s} = [];
                continue
            end

            % 1 -> 0 deg, 11 -> 72 deg, etc.
            phi = mod((phiIdx - 0.5) * (2*pi/50), 2*pi);
    end
    % figure; 
    % subplot 211
    % histogram(phi, [0:pi/5:2*pi])
    % subplot 212
    % histogram(phiIdx, [0:50/10:50])
    phiCell{s}  = phi(:);
    nBreaths(s) = numel(phi);

    m = mean(exp(1i*phi));
    Rvec(s)  = abs(m);
    muDeg(s) = mod(angle(m) * 180/pi, 360);
end

% ---------- histogram counts ----------
edges = linspace(0, 2*pi, nBins+1);
binCtrDeg = mod((edges(1:end-1) + diff(edges)/2) * 180/pi, 360);

countMat = zeros(nSub, nBins);
for s = 1:nSub
    if ~isempty(phiCell{s})
        countMat(s,:) = histcounts(phiCell{s}, edges);
    end
end

totCounts = sum(countMat, 1);
maxR = max(totCounts);
if maxR <= 0
    maxR = 1;
end

% ---------- theme colors ----------
circleCol = 0.55*labCol + 0.45*bgCol;
tickCol   = labCol;

% ---------- figure / axes ----------
hFig = figure( ...
    'Color', bgCol, ...
    'Position', [50 50 1100 850], ...
    'InvertHardcopy', 'off');

ax = axes('Parent', hFig);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');

ax.Color      = bgCol;
ax.LineWidth  = 2.2;
ax.FontSize   = 16;
ax.FontWeight = 'bold';
ax.FontName   = 'Dotum';
ax.XColor     = labCol;
ax.YColor     = labCol;
box(ax, 'off');

% ---------- reference circles ----------
tt = linspace(0, 2*pi, 500);

refLevels = unique(round([0.25 0.5 0.75 1.0] * maxR));
refLevels = refLevels(refLevels > 0);

for i = 1:numel(refLevels)
    rr = refLevels(i);
    plot(ax, rr*cos(tt), rr*sin(tt), '-', ...
        'Color', circleCol, ...
        'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
end

% label reference circles
refAng = 35 * pi/180;
for i = 1:numel(refLevels)
    rr = refLevels(i);
    [xr, yr] = pol2cart(refAng, rr);
    text(ax, xr, yr, sprintf('%d', rr), ...
        'Color', tickCol, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom');
end

% ---------- stacked circular histogram ----------
binWidth = edges(2) - edges(1);
angGap   = min(binWidth * 0.10, 3*pi/180);

for b = 1:nBins
    th1 = edges(b)   + angGap/2;
    th2 = edges(b+1) - angGap/2;

    if th2 <= th1
        continue
    end

    cumR = 0;
    for s = 1:nSub
        c = countMat(s,b);
        if c <= 0
            continue
        end

        r1 = cumR;
        r2 = cumR + c;
        th = linspace(th1, th2, 25);

        [xo, yo] = pol2cart(th, r2*ones(size(th)));
        [xi, yi] = pol2cart(fliplr(th), r1*ones(size(th)));

        p = patch(ax, [xo xi], [yo yi], cols(s,:), ...
            'EdgeColor', cols(s,:), ...
            'LineWidth', 1.0, ...
            'FaceAlpha', 0.95);
        p.HandleVisibility = 'off';

        cumR = r2;
    end
end

% ---------- outer boundary ----------
plot(ax, maxR*cos(tt), maxR*sin(tt), '-', ...
    'Color', circleCol, ...
    'LineWidth', 2.5, ...
    'HandleVisibility', 'off');

% ---------- circular tick marks / labels ----------
rTck1 = maxR * 1.02;
rTck2 = maxR * 1.10;
rLab  = maxR * 1.22;

for k = 1:numel(tickDeg)
    ang = tickDeg(k) * pi/180;

    [x1,y1] = pol2cart(ang, rTck1);
    [x2,y2] = pol2cart(ang, rTck2);
    plot(ax, [x1 x2], [y1 y2], '-', ...
        'Color', tickCol, ...
        'LineWidth', 4, ...
        'HandleVisibility', 'off');

    [xl,yl] = pol2cart(ang, rLab);

    degNow = mod(tickDeg(k), 360);
    if degNow == 0
        ha = 'left';  va = 'middle';
    elseif degNow == 180
        ha = 'right'; va = 'middle';
    elseif degNow > 0 && degNow < 180
        ha = 'center'; va = 'bottom';
    else
        ha = 'center'; va = 'top';
    end

    text(ax, xl, yl, tickText{k}, ...
        'HorizontalAlignment', ha, ...
        'VerticalAlignment', va, ...
        'Interpreter', 'none', ...
        'FontSize', 16, ...
        'FontWeight', 'bold', ...
        'FontName', 'Dotum', ...
        'Color', tickCol);
end

% ---------- legend ----------
hLeg = gobjects(nSub,1);
for s = 1:nSub
    hLeg(s) = plot(ax, nan, nan, '-', ...
        'Color', cols(s,:), ...
        'LineWidth', 8);
end

% if nSub > 0
%     lgd = legend(ax, hLeg, cellstr(uName), ...
%         'Interpreter', 'none', ...
%         'Box', 'off', ...
%         'Location', 'southeastoutside');
% 
%     lgd.TextColor = labCol;
%     lgd.Color     = bgCol;
%     lgd.FontName  = 'Dotum';
%     lgd.FontSize  = 13;
% end

% ---------- title ----------
% title(ax, sprintf('Breath-level stacked circular histogram (%s)', srcLabel), ...
%     'Color', labCol, ...
%     'FontSize', 18, ...
%     'FontWeight', 'bold', ...
%     'FontName', 'Dotum');

% ---------- limits ----------
lim = maxR * 1.38;
xlim(ax, [-lim lim]);
ylim(ax, [-lim lim]);

% ---------- subjects x summary table output ----------
repInfo = allSubIDs(repRow, 1:8);
subTab  = cell2table(repInfo, 'VariableNames', ...
    {'sessID','task','comboID','unusedChan','SessNum2','TYPE','GlobalID','subID2'});

subTab.nBreaths = nBreaths;
subTab.R        = Rvec;
subTab.muDeg    = muDeg;

end