function fig_final()
%FIG_FINAL  The three finalised grant figures.
%   Sized and weighted for heavy reduction in a grant document: large type,
%   thick lines, thick axes, no patient initials on the page.

P = ohrv_config();
outDir = P.work; figDir = P.figs;
[cmap, subs] = ohrv_colors();

FS   = 22;   % axis label / tick font
LW   = 4.0;  % data line width
AXLW = 2.6;  % axis line width
MS   = 13;   % marker size

S = readtable(fullfile(outDir,'report_sessions.csv'),'TextType','string');
D = readtable(fullfile(outDir,'report_changes.csv'),'TextType','string');
Q = D(D.from==1 & D.to==2,:);

%% ---------------- FIGURE 1: olfactory recovery ----------------
f = figure('Position',[60 60 760 640],'Color','w'); hold on
for k = 1:numel(subs)
    m = sortrows(S(S.subj==subs(k),:),'sess');
    y = m.olf; g = isfinite(y);
    if sum(g) < 2, continue; end
    plot(m.sess(g), y(g), '-o', 'Color', cmap(k,:), 'MarkerFaceColor', cmap(k,:), ...
        'MarkerEdgeColor','w', 'LineWidth', LW, 'MarkerSize', MS+2);
end
yline(0, ':', 'Color', [.55 .58 .57], 'LineWidth', 2);
set(gca,'XTick',1:3,'FontSize',FS,'LineWidth',AXLW,'Box','off','TickDir','out', ...
    'XColor',[.12 .14 .13],'YColor',[.12 .14 .13]);
xlim([0.85 3.15]);
xlabel('Session','FontSize',FS+2,'FontWeight','bold');
ylabel('Olfactory Function (z)','FontSize',FS+2,'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'FIG1_olfactory_recovery.png'), 'Resolution', 300);

%% ---------------- FIGURE 2: vagal outflow vs olfaction ----------------
f = figure('Position',[60 60 760 640],'Color','w'); hold on
x = Q.dOlf; y = Q.d_b_vol;
p = polyfit(x,y,1); xx = linspace(-0.85,1.75,20);
plot(xx, polyval(p,xx), '-', 'Color', [.42 .48 .46], 'LineWidth', LW-0.7);
for k = 1:numel(subs)
    i = find(Q.subj == subs(k));
    if isempty(i), continue; end
    plot(x(i), y(i), 'o', 'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor','w', ...
        'MarkerSize', MS+9, 'LineWidth', 2);
end
yline(0,':','Color',[.55 .58 .57],'LineWidth',2);
xline(0,':','Color',[.55 .58 .57],'LineWidth',2);
set(gca,'FontSize',FS,'LineWidth',AXLW,'Box','off','TickDir','out', ...
    'XColor',[.12 .14 .13],'YColor',[.12 .14 .13]);
xlabel('Change in Olfactory Function','FontSize',FS+2,'FontWeight','bold');
ylabel('Change in Vagal Outflow','FontSize',FS+2,'FontWeight','bold');
xlim([-0.85 1.8]); ylim([-0.30 0.30]);
text(-0.78, 0.265, sprintf('n = %d', height(Q)), 'FontSize',FS,'FontWeight','bold');
text(-0.78, 0.205, sprintf('r = %.2f', corr(x,y)), 'FontSize',FS,'FontWeight','bold');
text(-0.78, 0.145, sprintf('\\rho = %.2f', corr(x,y,'Type','Spearman')), ...
    'FontSize',FS,'FontWeight','bold');
exportgraphics(f, fullfile(figDir,'FIG2_vagal_vs_olfaction.png'), 'Resolution', 300);

%% ---------------- FIGURE 3: JH, depth effect after removing length ----------------
% Fit RSA ~ breath length within each session, then plot the residual RSA
% (ms) against breath depth. The slope of that residual relation is the
% depth effect with the duration effect already removed.
sess = {'250818_Dupi_NMH_JH_1','250818_Dupi_NMH_JH_2'};
scol = [0.55 0.62 0.68; 0.05 0.40 0.36];      % session 1 grey-blue, session 2 deep teal
f = figure('Position',[60 60 780 640],'Color','w'); hold on
h = gobjects(1,2); slp = nan(1,2); allV = [];
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

    b   = [ones(height(T),1) T.len] \ T.RR_max_min;      % RSA ~ breath length
    res = (T.RR_max_min - [ones(height(T),1) T.len]*b) * 1000;   % residual, ms
    v   = T.inhVol;

    % raw breaths, faint - context only
    scatter(v, res, 26, scol(s,:), 'filled', 'MarkerFaceAlpha', 0.16, 'MarkerEdgeColor','none');

    % fit + 95% confidence band on the mean
    md = fitlm(v, res); slp(s) = md.Coefficients.Estimate(2);
    vv = linspace(prctile(v,1), prctile(v,99), 100)';
    [yh, ci] = predict(md, vv);
    fill([vv; flipud(vv)], [ci(:,1); flipud(ci(:,2))], scol(s,:), ...
        'FaceAlpha', 0.28, 'EdgeColor','none');
    h(s) = plot(vv, yh, '-', 'Color', scol(s,:), 'LineWidth', LW+2.5);

    % depth-binned means +/- SEM - collapses the scatter so the slopes read
    nb = 6; edges = prctile(v, linspace(0,100,nb+1)); edges(end) = edges(end)+eps;
    for q = 1:nb
        m = v >= edges(q) & v < edges(q+1);
        if sum(m) < 8, continue; end
        mx = mean(v(m)); my = mean(res(m)); se = std(res(m))/sqrt(sum(m));
        plot([mx mx], [my-se my+se], '-', 'Color', scol(s,:), 'LineWidth', 3.4);
        plot(mx, my, 'o', 'MarkerFaceColor', scol(s,:), 'MarkerEdgeColor','w', ...
            'MarkerSize', 17, 'LineWidth', 2.2);
    end
    allV = [allV; v]; %#ok<AGROW>
    cb = coefCI(md);
    fprintf('%s: n=%d  depth slope = %+.3f ms/unit  (95%% CI %+.3f to %+.3f)\n', ...
        sess{s}, height(T), slp(s), cb(2,1), cb(2,2));
end
fprintf('slope difference (S2 - S1) = %+.3f ms per unit volume\n', slp(2)-slp(1));
yline(0,':','Color',[.55 .58 .57],'LineWidth',2);
% Clip to the central 98% of breath depths. A handful of very deep breaths
% otherwise stretch the axis and compress the bulk of the data into the left
% third of the panel; the fits are unchanged, only the visible range.
xlim([0 prctile(allV, 99)]);
ylim([-260 380]);
set(gca,'FontSize',FS,'LineWidth',AXLW,'Box','off','TickDir','out', ...
    'XColor',[.12 .14 .13],'YColor',[.12 .14 .13]);
xlabel('Breath Depth (inhaled volume, a.u.)','FontSize',FS+2,'FontWeight','bold');
ylabel('respHRV Residual (ms)','FontSize',FS+2,'FontWeight','bold');
lg = legend(h, {'Session 1','Session 2'}, 'Location','northwest','FontSize',FS-1);
lg.Box = 'off';
exportgraphics(f, fullfile(figDir,'FIG3_JH_depth_residual.png'), 'Resolution', 300);

fprintf('\nFIG_FINAL DONE\n');
end
