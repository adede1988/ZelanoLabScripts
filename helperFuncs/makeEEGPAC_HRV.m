function fig = makeEEGPAC_HRV(chanDat, macChan, taskVec, conds)


x = 5

useVec = chanDat.use == 1; 

% ============================================================
% Panels:
%  (1) Polar: mean HRV vs Phi (binned), one curve per band
%  (2) Scatter: HRV (y) vs Gamma power (x)
%  (3) Scatter: PAC (y) vs Gamma power (x), 3 bands + fit lines + p summary
% ============================================================

fig = figure('Color','w', 'visible', true, 'position', [0,0,1200, 700]);
tlo = tiledlayout(fig, 2, 3, 'Padding','compact', 'TileSpacing','compact');

pacFrex = chanDat.pac.PACfrex;

freqBreaks = [2, 8; ...
              8, 14; ...
              14,30];
freqLabs = {'theta','alpha','beta'};

% colors: forest green, burnt orange, dark mauve
frexCols = [0.13 0.55 0.13;
            0.80 0.33 0.00;
            0.55 0.41 0.53];

% ---------- build per-breath vectors for each band ----------
HRV = double(chanDat.behDat.HRV_RSAamp(:));
gam = max(double(macChan.gammaEnv.gamEnv), [], 2);   % [nBreath x 1] (gamma "power" proxy)
use = logical(useVec(:));

nB = numel(use);
nBand = 3;

PAC = nan(nB, nBand);
PHI = nan(nB, nBand);

for b = 1:nBand
     fidx = (pacFrex > freqBreaks(b,1)) & (pacFrex <= freqBreaks(b,2));
    if ~any(fidx), continue; end

    m = use;                                % keep your masking style
    idxBreath = find(m);                    % absolute breath indices
    fList = find(fidx);                     % absolute freq indices in pacFrex

    % --- PAC strength per breath (max over phase and freqs) ---
    Pz = chanDat.pac.pac(m, 1:50, fidx, 10);  % [nUse x 50 x nFband]
    curpac = squeeze(median(max(Pz, [], 2, 'omitnan'), 3, 'omitnan')); % [nUse x 1]
    PAC(idxBreath,b) = double(curpac);


       % --- find (phase,freq) of peak PAC per breath (in the z slice) ---
    tmp = squeeze(max(Pz, [], 3, 'omitnan'));                 % [nUse x 50]
    [~, maxi] = max(tmp, [], 2, 'omitnan');                   % [nUse x 1], phase bin

    tmp = squeeze(max(Pz, [], 2, 'omitnan'));                 % [nUse x nFband]
    [~, maxf] = max(tmp, [], 2, 'omitnan');                   % [nUse x 1], index within band

    % --- YOUR WAY of pulling per-breath phase preference at that peak ---
    freqAbs = fList(maxf);                                    % map band-index -> absolute freq index
    phi = arrayfun(@(x,y,z) chanDat.pac.pac(z, x, y, 5), maxi, freqAbs, idxBreath);
    phi = double(phi(:));

    % make sure phi is radians in [-pi,pi]
    if ~isreal(phi), phi = angle(phi); end
    if max(abs(phi),[],'omitnan') > (2*pi + 0.25), phi = deg2rad(phi); end
    phi = mod(phi + pi, 2*pi) - pi;

    PHI(idxBreath,b) = phi;
end

% common mask
mBase = use & isfinite(HRV) & isfinite(gam);

% ============================================================
% (1) Polar plot: mean HRV vs Phi, one curve per band
% ============================================================
ax1_cart = nexttile(tlo, 1); pos = ax1_cart.Position; delete(ax1_cart);
ax1 = polaraxes('Position', pos); hold(ax1,'on');

nBins = 24;
edges = linspace(-pi, pi, nBins+1);
cent  = edges(1:end-1) + diff(edges)/2;

for b = 1:nBand
    mb = mBase & isfinite(PHI(:,b));
    if nnz(mb) < 10, continue; end

    phi = PHI(mb,b);
    hrv = HRV(mb);

    % bin means
    bin = discretize(phi, edges);
    r = nan(1,nBins);
    for k = 1:nBins
        r(k) = mean(hrv(bin==k), 'omitnan');
    end

    th = [cent, cent(1)];      % close curve
    rr = [r,    r(1)];

    polarplot(ax1, th, rr, 'LineWidth', 2, 'Color', frexCols(b,:));
end

title(ax1, 'Mean HRV vs \phi (binned)', 'FontWeight','bold');
ax1.LineWidth = 1.5;
% (polar axes don't support grid alpha the same way; this keeps it clean)
hold(ax1,'off');
legend(ax1, freqLabs, 'Location','southoutside');
% ============================================================
% (2) Scatter: HRV vs Gamma power
% ============================================================
ax2 = nexttile(tlo, 2); cla(ax2); hold(ax2,'on');

scatter(ax2, gam(mBase), HRV(mBase), 10, [0 0 0], 'filled', ...
    'MarkerFaceAlpha', 0.55, 'MarkerEdgeAlpha', 0.55);

xlabel(ax2, 'Gamma power (per breath)');
ylabel(ax2, 'HRV\_RSAamp');
title(ax2, 'HRV vs OB \gamma power', 'FontWeight','bold');

grid(ax2,'on'); ax2.GridAlpha = 0.15;
ax2.LineWidth = 1.5; ax2.FontSize = 11;
box(ax2,'off');
hold(ax2,'off');

% ============================================================
% (3) Scatter: PAC vs Gamma power (3 bands) + fits + p summary
% ============================================================
ax3 = nexttile(tlo, 3); cla(ax3); hold(ax3,'on');

pVals = nan(1,nBand);

for b = 1:nBand
    mb = mBase & isfinite(PAC(:,b));
    if nnz(mb) < 10, continue; end

    x = gam(mb);
    y = PAC(mb,b);

    h = scatter(ax3, x, y, 12, frexCols(b,:), 'filled');
    h.MarkerFaceAlpha = 0.65;
    h.MarkerEdgeAlpha = 0.65;

    mdl = fitlm(x, y);
    pVals(b) = mdl.Coefficients.pValue(2);

    xx = linspace(min(x), max(x), 200);
    yy = mdl.Coefficients.Estimate(1) + mdl.Coefficients.Estimate(2)*xx;
    plot(ax3, xx, yy, '-', 'Color', frexCols(b,:), 'LineWidth', 2, 'HandleVisibility','off');
end

% p summary text: p = ##, ##, ##
pStr = strings(1,nBand);
for b = 1:nBand
    if isfinite(pVals(b)), pStr(b) = sprintf('%.3g', pVals(b));
    else,                 pStr(b) = 'NA';
    end
end
text(ax3, 0.02, 0.98, sprintf('p = %s, %s, %s', pStr(1), pStr(2), pStr(3)), ...
    'Units','normalized', 'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'FontWeight','bold', 'BackgroundColor','w', 'Margin', 4);

xlabel(ax3, 'Gamma power (per breath)');
ylabel(ax3, 'PAC strength (z)');
title(ax3, 'PAC vs OB \gamma power', 'FontWeight','bold');

grid(ax3,'on'); ax3.GridAlpha = 0.15;
ax3.LineWidth = 1.5; ax3.FontSize = 11;
box(ax3,'off');
hold(ax3,'off');

% (optional) add a small legend just for band colors
% 

% ============================================================
% Drop-in: 3 panels comparing PAC strength across bands
%   theta vs alpha, alpha vs beta, theta vs beta
% Dots colored by HRV (uses PAC, HRV, mBase already computed)
% ============================================================

pairs = [1 2; 2 3; 1 3];
pairLabs = {'\theta vs \alpha','\alpha vs \beta','\theta vs \beta'};

% common mask: need PAC in both bands + HRV finite
mHRV = mBase & isfinite(HRV);

for k = 1:3
    i = pairs(k,1); j = pairs(k,2);
    ax = nexttile(tlo, k+3);  % assumes you already have a tiledlayout active

    m = mHRV & isfinite(PAC(:,i)) & isfinite(PAC(:,j));
    if nnz(m) < 10
        title(ax, sprintf('PAC %s (n=%d)', pairLabs{k}, nnz(m)));
        box(ax,'off');
        continue
    end

    scatter(ax, PAC(m,i), PAC(m,j), 20, HRV(m), 'filled');
    colormap(ax, parula);
    cb = colorbar(ax);
    cb.Label.String = 'HRV\_RSAamp';
    clim(ax, prctile(HRV(m), [10 90])); 

    xlabel(ax, sprintf('PAC %s (z)', freqLabs{i}));
    ylabel(ax, sprintf('PAC %s (z)', freqLabs{j}));
    title(ax, sprintf('PAC: %s (colored by HRV)', pairLabs{k}));

    grid(ax,'on'); ax.GridAlpha = 0.15;
    ax.LineWidth = 1.5; ax.FontSize = 11;
    box(ax,'off');

    % optional: equal axis scaling for easier comparison
    axis(ax,'square');
end









end