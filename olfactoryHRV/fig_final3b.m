function fig_final3b(subjCode, sessNums, nBins, heightMult, fontScale)
%FIG_FINAL3B  Notched box plots of respHRV residual by breath-depth quantile.
%   fig_final3b()              -> JH, sessions 1 and 2, tertiles
%   fig_final3b('PC')          -> PC, sessions 1 and 2, tertiles
%   fig_final3b('JH',[1 2],5)  -> JH, quintiles
%   fig_final3b('AB',[1 2])    -> any patient / session pair present in work/
%
%   respHRV ~ breath length is fitted within each session and the residual (ms)
%   is summarised in NBINS equal-count breath-depth bins (default 3 = tertiles).
%   The output name gets a _q<N> suffix for any NBINS other than 3, so trying a
%   different binning does not overwrite the tertile figure.
%
%   Quantile edges are computed from BOTH sessions POOLED within a patient, then
%   applied to each session. That matters: per-session quantiles would put
%   different absolute breath depths into the same bin and the paired comparison
%   would stop being like-for-like. Pooling also gives roughly equal counts per
%   box, which the earlier fixed-width bins did not (PC had n=5 in the shallow
%   bin). Consequence: edges differ BETWEEN patients, so one patient's x axis is
%   not directly comparable to another's - the console prints the actual ranges.
%
%   Each bin holds a touching pair of boxes, earlier session left and later
%   session right, with the x tick centred on the pair. Boxes are notched: full
%   width at the quartiles and indented to the median, where the notch spans the
%   95% CI of the median (median +/- 1.57*IQR/sqrt(n)), so non-overlapping
%   notches indicate a reliable median difference. No whiskers: individual
%   breaths are drawn as a jittered cloud on top of a translucent box, so the raw
%   spread is read from the points themselves rather than from a percentile bar.
%
%   Drawn by hand rather than with boxplot() to control pairing, notch geometry
%   and draw order.

if nargin < 1 || isempty(subjCode),   subjCode   = 'JH'; end
if nargin < 2 || isempty(sessNums),   sessNums   = [1 2]; end
if nargin < 3 || isempty(nBins),      nBins      = 3;     end
if nargin < 4 || isempty(heightMult), heightMult = 1;     end
if nargin < 5 || isempty(fontScale),  fontScale  = 1;     end
subjCode = char(subjCode);

P = ohrv_config();
outDir = P.work; figDir = P.figs;

% Type and line weights scale with the CANVAS WIDTH, not with a fixed point
% size. When a figure is shrunk to fit a column, the rendered font size is
% FS * (columnWidth / figureWidth) - so a wider canvas at a fixed FS produces
% smaller text on the page. Scaling FS with width keeps text the same size on
% the page across the 3-bin and 5-bin versions; FONTSCALE bumps it further.
FIGW = max(900, 250*nBins+150);
K    = (FIGW/900) * fontScale;
FS   = round(22 * K);   % axis font
AXLW = 2.6 * K;         % axis line width
BOXLW = 2.4 * K;        % box edge
MEDLW = 3.6 * K;        % median bar
MSZ   = 22  * K;        % scatter marker area
CENT = 1:nBins;
BW   = 0.44;   % width of a single box in axis units; the pair touches at CENT
IND  = 0.34;   % horizontal indent at the median, as a fraction of BW

switch nBins
    case 3, LAB = {'Bottom third','Middle third','Top third'};
            XL  = 'Breath Depth';
    case 4, LAB = {'1st','2nd','3rd','4th'};
            XL  = 'Breath Depth (quartile)';
    case 5, LAB = {'1st','2nd','3rd','4th','5th'};
            XL  = 'Breath Depth (quintile)';
    otherwise
            LAB = arrayfun(@(k) sprintf('%d', k), 1:nBins, 'uni', 0);
            XL  = sprintf('Breath Depth (%d-quantile)', nBins);
end

% locate this patient's slim extracts
sess = cell(1,2);
for s = 1:2
    d = dir(fullfile(outDir, sprintf('*_%s_%d_slim.mat', subjCode, sessNums(s))));
    if isempty(d)
        error('fig_final3b:noData', 'no slim extract for %s session %d in %s', ...
              subjCode, sessNums(s), outDir);
    end
    sess{s} = erase(d(1).name, '_slim.mat');
end
fprintf('%s: %s  vs  %s\n', subjCode, sess{1}, sess{2});

% ---- pass 1: load both sessions, fit out breath length, keep residuals ----
S = struct('v',{[],[]}, 'res',{[],[]}, 'slope',{NaN,NaN});
for s = 1:2
    L = load(fullfile(outDir,[sess{s} '_slim.mat']));
    T = L.T; fs = L.fs; nSamp = L.nSamp; clear L
    keep = T.goodBreath==1 & isfinite(T.RR_max_min) & T.RR_max_min>0 ...
         & T.len>=1.5 & T.len<=15 ...
         & isfinite(T.inhDur) & T.inhDur>0 & T.inhDur<T.len ...
         & isfinite(T.inhVol) & T.inhVol>0 ...
         & isfinite(T.finalOnset) & T.finalOnset>=1 ...
         & (T.finalOnset+round(T.len*fs))<=nSamp & ~strcmpi(T.noseMouth,"mouth");
    T = T(keep,:);

    b = [ones(height(T),1) T.len] \ T.RR_max_min;                       % respHRV ~ breath length
    S(s).res   = (T.RR_max_min - [ones(height(T),1) T.len]*b) * 1000;   % residual, ms
    S(s).v     = T.inhVol;
    pf         = polyfit(S(s).v, S(s).res, 1);
    S(s).slope = pf(1);
end

% ---- quantile edges from both sessions pooled ----
allV  = [S(1).v; S(2).v];
cuts  = linspace(0, 100, nBins+1);
EDGES = [min(allV) prctile(allV, cuts(2:end-1)) max(allV)];
EDGES(end) = EDGES(end) + eps(EDGES(end));
fprintf('  depth edges (both sessions pooled):%s  a.u.\n', sprintf(' %.0f |', EDGES));

% ---- pass 2: draw ----
scol = [0.55 0.62 0.68; 0.05 0.40 0.36];
ecol = [0.30 0.38 0.45; 0.02 0.24 0.21];   % darker edge / median

rng(3);   % reproducible jitter
% HEIGHTMULT stretches the axes vertically without touching the width, so the
% same data occupies more vertical space and small median differences separate.
f = figure('Position',[60 60 FIGW round(660*heightMult)], ...
           'Color','w'); hold on
h = gobjects(1,2);

for s = 1:2
    v = S(s).v; res = S(s).res;
    for q = 1:nBins
        m = v >= EDGES(q) & v < EDGES(q+1);
        n = sum(m);
        if n < 5
            fprintf('  S%d  %-13s n=%3d  (too few, box omitted)\n', s, LAB{q}, n);
            continue
        end
        r = res(m);
        q1 = prctile(r,25); q2 = median(r); q3 = prctile(r,75);
        iqr = q3 - q1;
        nh  = 1.57*iqr/sqrt(n);                       % notch half-height
        nl  = max(q2-nh, q1); nu = min(q2+nh, q3);    % keep the notch inside the box

        % box column: earlier session left of the tick, later session right
        if s == 1, xl = CENT(q)-BW; xr = CENT(q); else, xl = CENT(q); xr = CENT(q)+BW; end
        d = BW*IND;

        % notched box: wide at the quartiles, indented to the median.
        % Drawn FIRST and kept translucent so the individual breaths sit on
        % top of it rather than under it.
        px = [xl xl xl+d xl xl xr xr xr-d xr xr];
        py = [q1 nl q2    nu q3 q3 nu q2    nl q1];
        patch(px, py, scol(s,:), 'FaceAlpha', 0.30, 'EdgeColor', ecol(s,:), ...
              'LineWidth', BOXLW);
        plot([xl+d xr-d],[q2 q2],'-','Color',ecol(s,:),'LineWidth',MEDLW);

        % jittered individual breaths, drawn on top of the box
        jx = (xl+xr)/2 + (rand(n,1)-0.5)*BW*0.80;
        scatter(jx, r, MSZ, scol(s,:), 'filled', 'MarkerFaceAlpha', 0.60, ...
                'MarkerEdgeColor','none');

        fprintf('  S%d  %-13s n=%3d  med=%+6.1f  IQR=%5.1f  notch=+/-%.1f\n', ...
            s, LAB{q}, n, q2, iqr, nh);
    end
    h(s) = patch(NaN, NaN, scol(s,:), 'FaceAlpha', 0.30, 'EdgeColor', ecol(s,:), ...
                 'LineWidth', BOXLW);   % legend proxy
end
fprintf('  depth slope: S%d %+.3f -> S%d %+.3f ms/unit\n', ...
    sessNums(1), S(1).slope, sessNums(2), S(2).slope);

yline(0,':','Color',[.55 .58 .57],'LineWidth',2);
set(gca,'FontSize',FS,'LineWidth',AXLW,'Box','off','TickDir','out', ...
    'XTick',CENT,'XTickLabel',LAB,'XTickLabelRotation',0, ...
    'XColor',[.12 .14 .13],'YColor',[.12 .14 .13]);
xlim([0.42 nBins+0.58]); ylim([-230 300]);
xlabel(XL,'FontSize',FS+2,'FontWeight','bold');
ylabel('respHRV Residual (ms)','FontSize',FS+2,'FontWeight','bold');
lg = legend(h, arrayfun(@(k) sprintf('Session %d',k), sessNums, 'uni', 0), ...
            'Location','northwest','FontSize',FS-1);
lg.Box = 'off';

suffix = '';
if nBins ~= 3,      suffix = sprintf('%s_q%d', suffix, nBins);      end
if heightMult ~= 1, suffix = sprintf('%s_h%g', suffix, heightMult); end
fn = fullfile(figDir, sprintf('FIG3b_%s_depth_boxes%s.png', subjCode, suffix));
exportgraphics(f, fn, 'Resolution', 300);
fprintf('wrote %s\n', fn);
end
