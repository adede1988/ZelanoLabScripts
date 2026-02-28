%% Group mean powZ heatmaps (assumes allChan is already in workspace)
% Computes per-channel meanPowZ = mean(tf.powZ(use==1,:,:),1,'omitnan')
% then plots group-mean heatmaps for requested subsets.
%
% MINIMAL CHANGE: adds plot_diff() calls + helper to make Dupi sess2-sess1 difference maps
% without changing how your existing plot_group() frequency axis labeling works.

% ---------- per-channel mean(powZ over breaths) ----------
n = numel(allChan);
chanMeanPowZ = cell(n,1);
chanMeanB = cell(n,1);
tLen = nan(n,1);
fLen = nan(n,1);

for i = 1:n
    powZ = allChan(i).tf.powZ;                 % breaths x time x frex
    breathSeg = allChan(i).tf.breathSeg;
    use  = (allChan(i).use(:) == 1);           % breaths x 1 logical

    %quick and dirty selection for audio book only in breathing task
    try 
        use = use & cellfun(@(x)strcmp(x, 'focus'), allChan(i).behDat.task); 
        disp('here')
    catch
    end
    m = powZ(use,:,:);
    b = breathSeg(use,:);

    extremeVals = arrayfun(@(x) sum(m(x,:,:)>100, 'all'), 1:sum(use));
    m = m(extremeVals < 100, :, :);
    b = b(extremeVals < 100, :);

    m = squeeze(mean(m, 1, 'omitnan'));  % time x frex
    b = squeeze(mean(b, 1, 'omitnan'));  % time x 1 (or 1 x time)

    chanMeanPowZ{i} = m;
    chanMeanB{i} = b(:);
    tLen(i) = size(m,1);
    fLen(i) = size(m,2);
end

% ---------- group indices ----------
task    = string({allChan.task});
type    = string({allChan.type});
sessNum = [allChan.sessNum];

idxDupi1 = find(task=="breathingTask" & sessNum==1 & type=="Dupi");
idxDupi2 = find(task=="breathingTask" & sessNum==2 & type=="Dupi");
idxOBE   = find(task=="breathingTask" & type=="OBE");   % ignore sessNum

idxDupCue1 = find(task=="cueTask" & sessNum==1 & type=="Dupi");
idxDupCue2 = find(task=="cueTask" & sessNum==2 & type=="Dupi");
idxOBECue  = find(task=="cueTask" & type=="OBE");       % ignore sessNum

idxDupO151 = find(task=="O15" & sessNum==1 & type=="Dupi");
idxDupO152 = find(task=="O15" & sessNum==2 & type=="Dupi");
idxOBEO15  = find(task=="O15" & type=="OBE");           % ignore sessNum

idxDupThresh1 = find(task=="threshTask" & sessNum==1 & type=="Dupi");
idxDupThresh2 = find(task=="threshTask" & sessNum==2 & type=="Dupi");

% ---------- make + plot group means ----------
plot_group(idxDupi1, "focusBreathing | sessNum=1 | type=Dupi");
plot_group(idxDupi2, "focusBreathing | sessNum=2 | type=Dupi");
plot_group(idxOBE,   "focusBreathing | type=Control");
% 
% plot_group(idxDupCue1, "cueTask | sessNum=1 | type=Dupi");
% plot_group(idxDupCue2, "cueTask | sessNum=2 | type=Dupi");
% plot_group(idxOBECue,  "cueTask | type=Control");
% 
% plot_group(idxDupO151, "O15 | sessNum=1 | type=Dupi");
% plot_group(idxDupO152, "O15 | sessNum=2 | type=Dupi");
% plot_group(idxOBEO15,  "O15 | type=Control");


% plot_group(idxDupThresh1, "Thresh | sessNum=1 | type=Dupi");
% plot_group(idxDupThresh2, "Thresh | sessNum=2 | type=Dupi");

% ---------- MINIMAL ADD: Dupi sess2 - sess1 difference heatmaps ----------
% plot_diff(idxDupi1,   idxDupi2,   "breathingTask diff(sess2_minus_sess1) | type=Dupi");
% plot_diff(idxDupCue1, idxDupCue2, "cueTask diff(sess2_minus_sess1) | type=Dupi");
% plot_diff(idxDupO151, idxDupO152, "O15 diff(sess2_minus_sess1) | type=Dupi");
% plot_diff(idxDupThresh1, idxDupThresh2, "Thresh diff(sess2_minus_sess1) | type=Dupi");

%% ---------- local helper ----------
function plot_group(idx, titleStr)
    if isempty(idx)
        warning("No entries found for: %s", titleStr);
        return
    end

    % Pull shared vars from base workspace (quick script style)
    allChan = evalin('base','allChan');
    chanMeanB = evalin('base','chanMeanB');
    chanMeanPowZ = evalin('base','chanMeanPowZ');
    tLen = evalin('base','tLen');
    fLen = evalin('base','fLen');

    % crop to common size (min across selected)
    tMin = min(tLen(idx));
    fMin = min(fLen(idx));

    M = nan(tMin, fMin, numel(idx));
    B = nan(tMin, numel(idx));
    for k = 1:numel(idx)
        tmp = chanMeanPowZ{idx(k)};
        M(:,:,k) = tmp(1:tMin, 1:fMin);
        tmp = chanMeanB{idx(k)};
        B(:,k) = tmp(1:tMin);
    end

    groupMean = mean(M, 3, 'omitnan'); % time x frex
    groupB    = mean(B, 2, 'omitnan');

    % frequency axis (best effort; assumes consistent frex ordering)
    frex = allChan(idx(1)).tf.frex(:);
    frex = frex(1:min(numel(frex), fMin));

    % ---------------- PLOT (unchanged style) ----------------
    figure('Color','w', 'position', [0,0,600,400]);

    % ---- limit freqs to >= 2 Hz ----
    fMask   = frex >= 2;
    frexUse = frex(fMask);

    % groupMean is [time x frex]
    gm = groupMean(:, fMask);                 % [tMin x nF>=2]

    % ---- heatmap (freq x time) ----
    imagesc(1:tMin, [], gm');            % y-axis = index into frexUse (YOUR STYLE)
    axis xy
    clim([-3 4])

    cb = colorbar;
    cb.Label.String = 'power (z-score)';

    % ---- axes styling (darker, bolder ticks) ----
    ax = gca;
    ax.Box       = 'off';
    ax.LineWidth = 1.5;
    ax.FontWeight= 'bold';
    ax.FontSize  = 12;
    ax.XColor    = [0 0 0];
    ax.YColor    = [0 0 0];
    ax.TickDir   = 'out';

    % ---- x-axis: labels only, no tick marks ----
    xticks([5 15 25 35 45])
    xticklabels({'inhale rise','inhale fall','exhale rise','exhale fall','pause'})
    ax.TickLength = [0 0];                    % removes tick marks but keeps labels
    xlabel('normalized time across sniff')

    % ---- y-axis (left) ----
    yyaxis left
    ylabel('Frequency (Hz)')
    yticks(20:20:180)
    yticklabels(round(frexUse(20:20:180)))

    title(sprintf('%s  (n=%d)', titleStr, numel(idx)), 'Interpreter','none')

    % ---- vertical dotted epoch boundaries ----
    hold on
    xline(10.5, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
    xline(20.5, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
    xline(30.5, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
    xline(40.5, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 2);

    % ---- respiration overlay on right axis ----
    yyaxis right
    if exist('groupB','var') && ~isempty(groupB)
        if numel(groupB) ~= tMin
            xOld = linspace(1, tMin, numel(groupB));
            groupB = interp1(xOld, groupB(:), 1:tMin, 'linear', 'extrap');
        end
        plot(1:tMin, groupB, 'k--', 'LineWidth', 3);   % thicker + dashed
    end
    ylabel('normalized respiration (inhale up)')
    ax.YAxis(2).Color = [0 0 0];
    ax.YAxis(2).LineWidth = 1.5;

    hold off

    % ---------------- SAVE (unchanged) ----------------
    saveDir = "R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs";
    if ~exist(saveDir,'dir'), mkdir(saveDir); end

    taskTok = regexp(titleStr, '^\s*([^|]+?)\s*(?:\||$)', 'tokens', 'once');
    if isempty(taskTok), taskTok = {char(titleStr)}; end
    taskName = string(strtrim(taskTok{1}));

    sessTok = regexp(titleStr, 'sessNum\s*=\s*(\d+)', 'tokens', 'once');
    if isempty(sessTok)
        sessPart = "sessAll";
    else
        sessPart = "sess" + string(sessTok{1});
    end

    typeTok = regexp(titleStr, 'type\s*=\s*([^|]+)', 'tokens', 'once');
    if isempty(typeTok)
        typePart = "typeNA";
    else
        typePart = "type" + string(strtrim(typeTok{1}));
    end

    safe = @(s) regexprep(string(s), '[^\w\-]+', '_');
    baseName = safe(taskName) + "_" + safe(sessPart) + "_" + safe(typePart) + "_powZ";

    outFile = fullfile(saveDir, baseName + ".jpg");
    exportgraphics(gcf, outFile, "Resolution", 600);
end


% ---------- MINIMAL ADD: difference plot helper (reuses your plotting style) ----------
function plot_diff(idxSess1, idxSess2, titleStr)
    if isempty(idxSess1) || isempty(idxSess2)
        warning("No entries found for diff: %s", titleStr);
        return
    end

    allChan = evalin('base','allChan');
    chanMeanB = evalin('base','chanMeanB');
    chanMeanPowZ = evalin('base','chanMeanPowZ');
    tLen = evalin('base','tLen');
    fLen = evalin('base','fLen');

    % crop to common size across BOTH groups
    tMin = min([min(tLen(idxSess1)) min(tLen(idxSess2))]);
    fMin = min([min(fLen(idxSess1)) min(fLen(idxSess2))]);

    % group 1 (sess1)
    M1 = nan(tMin, fMin, numel(idxSess1));
    B1 = nan(tMin, numel(idxSess1));
    for k = 1:numel(idxSess1)
        tmp = chanMeanPowZ{idxSess1(k)};
        M1(:,:,k) = tmp(1:tMin, 1:fMin);
        tmpB = chanMeanB{idxSess1(k)};
        B1(:,k) = tmpB(1:tMin);
    end
    g1 = mean(M1, 3, 'omitnan');
    b1 = mean(B1, 2, 'omitnan');

    % group 2 (sess2)
    M2 = nan(tMin, fMin, numel(idxSess2));
    B2 = nan(tMin, numel(idxSess2));
    for k = 1:numel(idxSess2)
        tmp = chanMeanPowZ{idxSess2(k)};
        M2(:,:,k) = tmp(1:tMin, 1:fMin);
        tmpB = chanMeanB{idxSess2(k)};
        B2(:,k) = tmpB(1:tMin);
    end
    g2 = mean(M2, 3, 'omitnan');
    b2 = mean(B2, 2, 'omitnan');


    % do stats: 
    alpha = 0.05;          % family-wise q threshold for FDR (BH)
    vartype = 'unequal';   % 'unequal' = Welch; set to 'equal' if you want pooled-variance t-test
    
    % ---- sanity checks ----
    [s1a,s1b,nSub1] = size(M1);
    [s2a,s2b,nSub2] = size(M2);
    
    
    nTests = s1a*s1b;
    
    % ---- reshape to [subjects x tests] so we can vectorize ----
    X = reshape(permute(M1,[3 1 2]), nSub1, nTests);   % [nSub1 x nTests]
    Y = reshape(permute(M2,[3 1 2]), nSub2, nTests);   % [nSub2 x nTests]
    
    % ---- vectorized two-sample t-tests across subjects (Dim=1) ----
    % NOTE: If you have NaNs that differ by pixel across subjects, see the NaN-safe loop below.
    [~, pVec, ~, stats] = ttest2(X, Y, 'Dim', 1, 'Vartype', vartype);
    
    tVec  = stats.tstat;            % 1 x nTests
    dfVec = stats.df;               % 1 x nTests (may vary for Welch)
    
    % ---- FDR (Benjamini–Hochberg) on all p-values ----
    qVec = nan(size(pVec));
    valid = isfinite(pVec);
    qVec(valid) = fdr_bh_q(pVec(valid));   % BH-adjusted p-values (aka q-values)
    
    % ---- reshape back to [50 x 300] maps ----
    pMap  = reshape(pVec,  [s1a s1b]);
    qMap  = reshape(qVec,  [s1a s1b]);
    tMap  = reshape(tVec,  [s1a s1b]);
    dfMap = reshape(dfVec, [s1a s1b]);
    
    sig_unc  = pMap < alpha;    % uncorrected
    sig_fdr  = qMap < alpha;    % FDR-corrected (BH)
    
    % ---- (optional) effect map: mean difference (M1 - M2) ----
    meanDiff = reshape(mean(X,1,'omitnan') - mean(Y,1,'omitnan'), [s1a s1b]);
    
    % Results you likely care about:
    %   tMap, pMap, qMap, sig_fdr, meanDiff, dfMap


    % difference: sess2 - sess1
    groupMean = g2 - g1;
    groupB    = mean([b1 b2], 2, 'omitnan');   % overlay respiration = average of the two

    % frequency axis (same approach as your plot_group)
    frex = allChan(idxSess1(1)).tf.frex(:);
    frex = frex(1:min(numel(frex), fMin));

    % ---------------- PLOT (copy of your existing plotting style) ----------------
    figure('Color','w', 'position', [0,0,600,400]);

    fMask   = frex >= 2;
    frexUse = frex(fMask);
    gm = groupMean(:, fMask);

    imagesc(1:tMin, [], gm');            % y-axis = index into frexUse (YOUR STYLE)
    axis xy
    clim([-2 2])                         % simple symmetric default for diffs

    cb = colorbar;
    cb.Label.String = 'power (z-score)';

    ax = gca;
    ax.Box       = 'off';
    ax.LineWidth = 1.5;
    ax.FontWeight= 'bold';
    ax.FontSize  = 12;
    ax.XColor    = [0 0 0];
    ax.YColor    = [0 0 0];
    ax.TickDir   = 'out';

    xticks([5 15 25 35 45])
    xticklabels({'inhale rise','inhale fall','exhale rise','exhale fall','pause'})
    ax.TickLength = [0 0];
    xlabel('normalized time across sniff')

    yyaxis left
    ylabel('Frequency (Hz)')
    yticks(20:20:180)
    yticklabels(round(frexUse(20:20:180)))

    title(sprintf('%s  (n2=%d, n1=%d)', titleStr, numel(idxSess2), numel(idxSess1)), 'Interpreter','none')

    hold on
    xline(10.5, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
    xline(20.5, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
    xline(30.5, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
    xline(40.5, ':', 'Color', [0.85 0.85 0.85], 'LineWidth', 2);

    yyaxis right
    if exist('groupB','var') && ~isempty(groupB)
        if numel(groupB) ~= tMin
            xOld = linspace(1, tMin, numel(groupB));
            groupB = interp1(xOld, groupB(:), 1:tMin, 'linear', 'extrap');
        end
        plot(1:tMin, groupB, 'k--', 'LineWidth', 3);
    end
    ylabel('normalized respiration (inhale up)')
    ax.YAxis(2).Color = [0 0 0];
    ax.YAxis(2).LineWidth = 1.5;

    hold off

    % ---------------- SAVE (same filename logic; diff lives in taskName) ----------------
    saveDir = "R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs";
    if ~exist(saveDir,'dir'), mkdir(saveDir); end

    taskTok = regexp(titleStr, '^\s*([^|]+?)\s*(?:\||$)', 'tokens', 'once');
    if isempty(taskTok), taskTok = {char(titleStr)}; end
    taskName = string(strtrim(taskTok{1}));

    sessTok = regexp(titleStr, 'sessNum\s*=\s*(\d+)', 'tokens', 'once');
    if isempty(sessTok)
        sessPart = "sessAll";
    else
        sessPart = "sess" + string(sessTok{1});
    end

    typeTok = regexp(titleStr, 'type\s*=\s*([^|]+)', 'tokens', 'once');
    if isempty(typeTok)
        typePart = "typeNA";
    else
        typePart = "type" + string(strtrim(typeTok{1}));
    end

    safe = @(s) regexprep(string(s), '[^\w\-]+', '_');
    baseName = safe(taskName) + "_" + safe(sessPart) + "_" + safe(typePart) + "_powZ";

    outFile = fullfile(saveDir, baseName + ".jpg");
    exportgraphics(gcf, outFile, "Resolution", 600);
end



function q = fdr_bh_q(p)
    p = p(:);
    m = numel(p);
    [ps, idx] = sort(p, 'ascend');
    qs = ps .* (m ./ (1:m)');              % raw BH
    qs = flipud(cummin(flipud(qs)));       % enforce monotonicity
    qs = min(qs, 1);
    q = nan(m,1);
    q(idx) = qs;
end