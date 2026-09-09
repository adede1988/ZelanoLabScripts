function fig_final3b(subjCode, sessNums)
%FIG_FINAL3B  Notched box plots of respHRV residual by breath depth, one patient.
%   fig_final3b()            -> JH, sessions 1 and 2 (the largest improver)
%   fig_final3b('PC')        -> PC, sessions 1 and 2 (the one decliner)
%   fig_final3b('AB',[1 2])  -> any patient / session pair present in work/
%
%   respHRV ~ breath length is fitted within each session and the residual (ms)
%   is summarised in three fixed 800-unit depth bins. Each bin holds a touching
%   pair of boxes, earlier session left and later session right, with the x tick
%   centred on the pair. Boxes are notched: full width at the quartiles and
%   indented to the median, where the notch spans the 95% CI of the median
%   (median +/- 1.57*IQR/sqrt(n)), so non-overlapping notches indicate a
%   reliable median difference. Individual breaths are overlaid as a jittered
%   cloud inside each box column.
%
%   Bins are held fixed across patients so figures are directly comparable.
%
%   Drawn by hand rather than with boxplot() to control pairing, notch geometry
%   and colours.

if nargin < 1 || isempty(subjCode), subjCode = 'JH'; end
if nargin < 2 || isempty(sessNums), sessNums = [1 2]; end
subjCode = char(subjCode);

P = ohrv_config();
outDir = P.work; figDir = P.figs;

FS   = 22;     % axis font
AXLW = 2.6;    % axis line width
EDGES = [0 800 1600 2400];
CENT  = EDGES(1:end-1) + diff(EDGES)/2;    % 400, 1200, 2000
BW    = 195;   % width of a single box; the pair spans 2*BW and touches at CENT
IND   = 0.34;  % horizontal indent at the median, as a fraction of BW

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

scol = [0.55 0.62 0.68; 0.05 0.40 0.36];
ecol = [0.30 0.38 0.45; 0.02 0.24 0.21];   % darker edge / median

rng(3);   % reproducible jitter
f = figure('Position',[60 60 900 660],'Color','w'); hold on
h = gobjects(1,2); slp = nan(1,2);

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

    b   = [ones(height(T),1) T.len] \ T.RR_max_min;              % respHRV ~ breath length
    res = (T.RR_max_min - [ones(height(T),1) T.len]*b) * 1000;   % residual, ms
    v   = T.inhVol;
    pf  = polyfit(v, res, 1); slp(s) = pf(1);

    for q = 1:numel(CENT)
        m = v >= EDGES(q) & v < EDGES(q+1);
        n = sum(m);
        if n < 5
            fprintf('  S%d  %4d-%4d  n=%3d  (too few, box omitted)\n', ...
                s, EDGES(q), EDGES(q+1), n);
            continue
        end
        r = res(m);
        q1 = prctile(r,25); q2 = median(r); q3 = prctile(r,75);
        w1 = prctile(r,10); w2 = prctile(r,90);
        iqr = q3 - q1;
        nh  = 1.57*iqr/sqrt(n);                       % notch half-height
        nl  = max(q2-nh, q1); nu = min(q2+nh, q3);    % keep the notch inside the box

        % box column: earlier session left of the tick, later session right
        if s == 1, xl = CENT(q)-BW; xr = CENT(q); else, xl = CENT(q); xr = CENT(q)+BW; end
        d = BW*IND;

        % jittered individual breaths, inside this column
        jx = (xl+xr)/2 + (rand(n,1)-0.5)*BW*0.80;
        scatter(jx, r, 15, scol(s,:), 'filled', 'MarkerFaceAlpha', 0.30, ...
                'MarkerEdgeColor','none');

        % whiskers
        xm = (xl+xr)/2;
        plot([xm xm],[w1 q1],'-','Color',ecol(s,:),'LineWidth',2.2);
        plot([xm xm],[q3 w2],'-','Color',ecol(s,:),'LineWidth',2.2);
        plot(xm+[-1 1]*BW*0.22,[w1 w1],'-','Color',ecol(s,:),'LineWidth',2.2);
        plot(xm+[-1 1]*BW*0.22,[w2 w2],'-','Color',ecol(s,:),'LineWidth',2.2);

        % notched box: wide at the quartiles, indented to the median
        px = [xl xl xl+d xl xl xr xr xr-d xr xr];
        py = [q1 nl q2    nu q3 q3 nu q2    nl q1];
        patch(px, py, scol(s,:), 'FaceAlpha', 0.80, 'EdgeColor', ecol(s,:), ...
              'LineWidth', 2.4);
        plot([xl+d xr-d],[q2 q2],'-','Color',ecol(s,:),'LineWidth',3.6);

        fprintf('  S%d  %4d-%4d  n=%3d  med=%+6.1f  IQR=%5.1f  notch=+/-%.1f\n', ...
            s, EDGES(q), EDGES(q+1), n, q2, iqr, nh);
    end
    h(s) = patch(NaN, NaN, scol(s,:), 'FaceAlpha', 0.80, 'EdgeColor', ecol(s,:), ...
                 'LineWidth', 2.4);   % legend proxy
end
fprintf('  depth slope: S%d %+.3f -> S%d %+.3f ms/unit\n', ...
    sessNums(1), slp(1), sessNums(2), slp(2));

yline(0,':','Color',[.55 .58 .57],'LineWidth',2);
set(gca,'FontSize',FS,'LineWidth',AXLW,'Box','off','TickDir','out', ...
    'XTick',CENT,'XTickLabel',{'0-800','800-1600','1600-2400'}, ...
    'XColor',[.12 .14 .13],'YColor',[.12 .14 .13]);
xlim([0 2400]); ylim([-230 300]);
xlabel('Breath Depth (inhaled volume, a.u.)','FontSize',FS+2,'FontWeight','bold');
ylabel('respHRV Residual (ms)','FontSize',FS+2,'FontWeight','bold');
lg = legend(h, arrayfun(@(k) sprintf('Session %d',k), sessNums, 'uni', 0), ...
            'Location','northwest','FontSize',FS-1);
lg.Box = 'off';
fn = fullfile(figDir, sprintf('FIG3b_%s_depth_boxes.png', subjCode));
exportgraphics(f, fn, 'Resolution', 300);
fprintf('wrote %s\n', fn);
end
