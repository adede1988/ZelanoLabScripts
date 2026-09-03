function rsa_report()
% Full downstream analysis: change scores for every slope and vagal measure,
% against olfactory change. Outcome-blind reliability gate first, then a
% max-statistic permutation over the retained family.

P      = ohrv_config();

sp     = P.code;
outDir = P.work;
figDir = P.figs;
rng(99); NPERT = 20000;

G = readtable(fullfile(outDir, 'session_gains.csv'), 'TextType', 'string');
O = readtable(fullfile(outDir, 'olfactory_scores.csv'), 'TextType', 'string');

%% ---------- olfactory scores ----------
O.thresh_mean = mean([O.thresh_low_cal, O.thresh_high_cal], 2, 'omitnan');
O.subj = extractBetween(O.sessID, "NMH_", "_" + digitsPattern);
O.sess = double(extractAfter(O.sessID, "NMH_" + lettersPattern + "_"));
O.subj(O.subj == "TPB") = "TB";
z = @(v) (v - mean(v, 'omitnan')) ./ std(v, 'omitnan');
O.olf = mean([z(O.cue_d), z(O.thresh_mean), z(O.O15_acc)], 2, 'omitnan');

M = innerjoin(G, O(:, {'subj','sess','cue_d','thresh_mean','O15_acc','olf'}), ...
              'Keys', {'subj','sess'});
M = sortrows(M, {'subj','sess'});
fprintf('=== %d sessions with gains + olfaction, %d subjects ===\n\n', ...
    height(M), numel(unique(M.subj)));
disp(M(:, {'sessID','subj','sess','n','nBlk','b_len','nullShift','nullIID','theta', ...
           'b_vol','nullVol','thetaVol','medRSA','HR','cue_d','O15_acc','olf'}));

%% ---------- change scores over adjacent intervals ----------
SLOPE = {'theta','thetaVol','b_len','b_vol','medRSA','HR'};
SE    = {'seTheta','seThetaVol','seTheta','seThetaVol','', ''};
D = table();
us = unique(M.subj);
for k = 1:numel(us)
    m = M(M.subj == us(k), :);
    for j = 1:height(m)-1
        row = table(us(k), m.sess(j), m.sess(j+1), ...
            m.olf(j+1)-m.olf(j), m.cue_d(j+1)-m.cue_d(j), m.O15_acc(j+1)-m.O15_acc(j), ...
            'VariableNames', {'subj','from','to','dOlf','dCue','dO15'});
        for q = 1:numel(SLOPE)
            row.(['d_' SLOPE{q}]) = m.(SLOPE{q})(j+1) - m.(SLOPE{q})(j);
            if ~isempty(SE{q})
                row.(['se_' SLOPE{q}]) = sqrt(m.(SE{q})(j)^2 + m.(SE{q})(j+1)^2);
            end
        end
        D = [D; row]; %#ok<AGROW>
    end
end

% merge the vagal panel deltas if present
pf = fullfile(outDir, 'panel_pairs.csv');
PAN = {};
if exist(pf, 'file')
    P = readtable(pf, 'TextType', 'string');
    P.subj = string(P.subj);
    D = outerjoin(D, P, 'Keys', {'subj','from','to'}, 'MergeKeys', true, 'Type', 'left');
    PAN = {'logRMSSD','logMASD','logMedPV','logHF','pNN20','evokedRSA','depthSlope','adjLogRSA'};
    PAN = PAN(ismember(cellfun(@(x) ['d_' x], PAN, 'uni', 0), D.Properties.VariableNames));
end
D = D(isfinite(D.dOlf), :);
writetable(D, fullfile(outDir, 'report_changes.csv'));

fprintf('\n=== change scores: %d intervals, %d subjects ===\n', height(D), numel(unique(D.subj)));
disp(D(:, [{'subj','from','to','dOlf','dCue','dO15'}, cellfun(@(x) ['d_' x], SLOPE, 'uni', 0)]));

%% ---------- primary contrast ----------
Q = D(D.from == 1 & D.to == 2, :);
fprintf('\nPRIMARY (session 1->2, matched interval): n = %d  [%s]\n', ...
    height(Q), strjoin(cellstr(Q.subj)', ', '));
fprintf('ALL adjacent intervals: n = %d\n', height(D));

%% ---------- reliability gate, outcome-blind ----------
CAND = [SLOPE(1:4), PAN];
fprintf('\n=== Gate 1: kappa = mean SE^2 / var(delta)  (outcome-blind) ===\n');
keepM = {};
for k = 1:numel(CAND)
    dcol = ['d_' CAND{k}]; scol = ['se_' CAND{k}];
    if ~ismember(dcol, Q.Properties.VariableNames), continue; end
    d = Q.(dcol);
    if ismember(scol, Q.Properties.VariableNames)
        se = Q.(scol); kap = mean(se.^2,'omitnan') / var(d,'omitnan');
    else
        se = nan(size(d)); kap = NaN;
    end
    ok = isfinite(kap) && kap <= 0.5;
    fprintf('  %-12s var=%9.4f  meanSE=%8.4f  kappa=%7.2f  %s\n', ...
        CAND{k}, var(d,'omitnan'), mean(se,'omitnan'), kap, tern(ok,'RETAIN','drop'));
    if ok, keepM{end+1} = CAND{k}; end %#ok<AGROW>
end
fprintf('\nretained: %s\n', tern(~isempty(keepM), strjoin(keepM, ', '), 'NONE'));

%% ---------- association with olfactory change ----------
for scope = ["primary", "all"]
    if scope == "primary", W = Q; else, W = D; end
    if height(W) < 4, continue; end
    fprintf('\n=== association with dOlf  [%s, n=%d] ===\n', scope, height(W));
    fprintf('%-12s %8s %8s %20s %10s\n', 'measure','rho','pearson','jackknife rho','P(rho>0)');
    rr = [];
    for k = 1:numel(keepM)
        dcol = ['d_' keepM{k}]; scol = ['se_' keepM{k}];
        d = W.(dcol); good = isfinite(d) & isfinite(W.dOlf);
        if sum(good) < 4, continue; end
        r = corr(d(good), W.dOlf(good), 'Type', 'Spearman');
        rp = corr(d(good), W.dOlf(good));
        jk = nan(sum(good),1); idx = find(good);
        for i = 1:numel(idx)
            m = good; m(idx(i)) = false;
            jk(i) = corr(d(m), W.dOlf(m), 'Type', 'Spearman');
        end
        if ismember(scol, W.Properties.VariableNames)
            se = W.(scol)(good); pr = nan(NPERT,1);
            for b = 1:NPERT
                pr(b) = corr(d(good) + se.*randn(sum(good),1), W.dOlf(good), 'Type','Spearman');
            end
            pgt = mean(pr > 0.01);
        else
            pgt = NaN;
        end
        fprintf('%-12s %+8.3f %+8.3f   [%+6.2f, %+6.2f]  %10.3f\n', ...
            keepM{k}, r, rp, min(jk), max(jk), pgt);
        rr(end+1) = r; %#ok<AGROW>
    end

    % max-statistic permutation, subject-level labels
    if ~isempty(rr)
        subs = unique(W.subj, 'stable');
        Pm = perms(1:numel(subs)); nP = size(Pm,1);
        Tobs = max(abs(rr)); Tn = nan(nP,1);
        for p = 1:nP
            map = containers.Map(cellstr(subs), cellstr(subs(Pm(p,:))));
            permOlf = W.dOlf;
            for i = 1:height(W)
                src = find(subs == string(map(char(W.subj(i)))), 1);
                sel = find(W.subj == subs(src), 1);
                permOlf(i) = W.dOlf(sel);
            end
            v = [];
            for k = 1:numel(keepM)
                d = W.(['d_' keepM{k}]); g = isfinite(d) & isfinite(permOlf);
                if sum(g) >= 4, v(end+1) = corr(d(g), permOlf(g), 'Type','Spearman'); end %#ok<AGROW>
            end
            if ~isempty(v), Tn(p) = max(abs(v)); end
        end
        pmax = mean(Tn >= Tobs - 1e-12, 'omitnan');
        fprintf('MAX |rho| = %.3f   subject-level permutation p = %.4f  (floor 1/%d = %.4f)\n', ...
            Tobs, pmax, nP, 1/nP);
    end

    % confounds
    fprintf('confounds:  rho(dOlf, dHR) = %+.3f', ...
        corr(W.d_HR, W.dOlf, 'Type', 'Spearman', 'Rows', 'complete'));
    if ismember('d_breathRate', W.Properties.VariableNames)
        fprintf('   rho(dOlf, dBreathRate) = %+.3f', ...
            corr(W.d_breathRate, W.dOlf, 'Type','Spearman','Rows','complete'));
    end
    fprintf('\n');
end

%% ---------- figures ----------
plotSet(Q, keepM, fullfile(figDir, 'report_primary.png'), 'session 1 to 2');
plotSet(D, keepM, fullfile(figDir, 'report_all.png'),     'all adjacent intervals');

f = figure('Position',[80 80 900 340],'Color','w');
subplot(1,2,1); hold on
errorbar(1:height(M), M.theta, M.seTheta, 'o', 'MarkerFaceColor',[.13 .47 .42], ...
    'Color',[.13 .47 .42], 'LineWidth',1.1, 'CapSize',5);
yline(0,':'); set(gca,'XTick',1:height(M),'XTickLabel',extractAfter(M.sessID,"NMH_"), ...
    'XTickLabelRotation',60,'TickLabelInterpreter','none','FontSize',7);
ylabel('\theta duration (corrected)'); title('Duration slope'); grid on; box on
subplot(1,2,2); hold on
errorbar(1:height(M), M.thetaVol, M.seThetaVol, 'o', 'MarkerFaceColor',[.64 .18 .28], ...
    'Color',[.64 .18 .28], 'LineWidth',1.1, 'CapSize',5);
yline(0,':'); set(gca,'XTick',1:height(M),'XTickLabel',extractAfter(M.sessID,"NMH_"), ...
    'XTickLabelRotation',60,'TickLabelInterpreter','none','FontSize',7);
ylabel('\theta depth (corrected)'); title('Depth slope'); grid on; box on
exportgraphics(f, fullfile(figDir,'report_sessions.png'), 'Resolution', 150);

writetable(M, fullfile(outDir, 'report_sessions.csv'));
fprintf('\nREPORT DONE\n');
end

function plotSet(W, keepM, fname, ttl)
if isempty(keepM) || height(W) < 3, return; end
nk = numel(keepM);
f = figure('Position',[60 60 min(300*nk,1500) 320],'Color','w');
for k = 1:nk
    subplot(1,nk,k); hold on
    d = W.(['d_' keepM{k}]);
    scol = ['se_' keepM{k}];
    if ismember(scol, W.Properties.VariableNames)
        errorbar(W.dOlf, d, W.(scol), 'o','MarkerFaceColor',[.13 .47 .42], ...
            'Color',[.13 .47 .42],'MarkerSize',7,'LineWidth',1.2,'CapSize',6);
    else
        plot(W.dOlf, d, 'o','MarkerFaceColor',[.13 .47 .42],'Color',[.13 .47 .42],'MarkerSize',7);
    end
    for i = 1:height(W)
        text(W.dOlf(i), d(i), "  " + W.subj(i) + string(W.from(i)) + ">" + string(W.to(i)), 'FontSize',7);
    end
    g = isfinite(d) & isfinite(W.dOlf);
    yline(0,':'); xline(0,':');
    title(sprintf('%s  \\rho=%+.2f', keepM{k}, corr(d(g),W.dOlf(g),'Type','Spearman')), ...
        'Interpreter','none','FontSize',9);
    xlabel('\Delta olfaction'); ylabel(['\Delta ' keepM{k}],'Interpreter','none');
    grid on; box on
end
sgtitle(ttl);
exportgraphics(f, fname, 'Resolution', 150);
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
