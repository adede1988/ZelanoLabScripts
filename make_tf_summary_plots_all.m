function make_tf_summary_plots_all()
% make_tf_summary_plots_all
% Batch-generate presentation-quality multi-panel TF/QC summary figures
% for all processed chanDat files in ONE folder (non-recursive).
%
% INPUTS (edit below):
%   inRoot  - folder containing processed chanDat .mat files
%   outRoot - folder to save jpg plots
%
% OUTPUT filename:
%   subID_sessNum_<chanLabel>_tfPlot.jpg
% where chanLabel = chanDat.labels{chanDat.chi} (sanitized for filesystem)
%
% Panels:
% 1) QC elimination proportions (4 reasons)
% 2) QC-passed breaths by condition (bars labeled task_noseMouth_warp)
% 3) Density of QC-passed breath lengths by condition
% 4) Mean respiration trace by condition (QC-passed)
% 5) Mean TF heatmap (powZ) log-freq axis + overlay mean breathSeg
% 6) Mean aperiodic_fit spectra per epoch (log-log)
% 7) Mean spectra_flat spectra per epoch (log-log)
% 8) Gamma peak frequency density per epoch (QC-passed)
%
% Notes:
% - Skips files that do not contain chanDat or required fields.
% - Uses ksdensity if available; otherwise histogram('pdf','stairs') fallback.

%% ---------------- USER SETTINGS ----------------
inRoot  = 'R:\Neurology\Zelano_Lab\Lab_Common\QuestMirror\CHANDAT_processed';
outRoot = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\SignalProcessing\tf_plots';

dpi = 300;                 % export resolution
visibleFigs = true;       % set true to watch plots generate

minBreathsPerCond = 1;     % minimum QC-passed breaths to plot a condition curve/trace

% TF heatmap color limits (set [] for auto)
tfCLim = [-5 10];

% Spectra scale (requested log-log)
spectraScale = 'xlinear_ylog';

% Keep raw respiration (requested)
useLowRspForPanel = false; %#ok<NASGU>  % left here in case you ever flip it

epochNames = {'inhale rise','inhale fall','exhale rise','exhale fall','pause'};
phaseEdges = [10.5 20.5 30.5 40.5];     % boundaries between 5 epochs (10 samples each)

%% ------------------------------------------------
if ~exist(outRoot,'dir'), mkdir(outRoot); end

% Non-recursive: only .mat in this folder
files = dir(fullfile(inRoot, '*.mat'));
fprintf('Found %d .mat files in %s\n', numel(files), inRoot);

for ii = 1:numel(files)
    fpath = fullfile(files(ii).folder, files(ii).name);
    fprintf('[%d/%d] %s\n', ii, numel(files), files(ii).name);

    try
        % quick check contains chanDat without loading everything
        w = whos('-file', fpath);
        if ~any(strcmp({w.name}, 'chanDat'))
            fprintf('  -> skipping (no chanDat)\n');
            continue
        end

        S = load(fpath, 'chanDat');
        chanDat = S.chanDat;

        % ---------------- metadata ----------------
        subID   = safeStr(getFieldOr(chanDat,'subID','UNK'));
        sessNum = getFieldOr(chanDat,'sessNum',1);
        if isempty(sessNum) || ~isfinite(sessNum), sessNum = 1; end
        sessNum = round(double(sessNum));

        chi = getFieldOr(chanDat,'chi',NaN);
        chi = double(chi);
        if ~isfinite(chi), chi = NaN; end

        chanLabel = getChanLabel(chanDat, chi);
        chanLabelClean = sanitizeForFilename(chanLabel);
        taskLab = chanDat.task; 
        typeLab = chanDat.type; 

        outName = sprintf('%s_%d_%s_%s_%s_tfPlot.jpg', subID, sessNum, chanLabelClean, taskLab, typeLab);
        outFile = fullfile(outRoot, outName);
if ~exist(outFile, 'file')
        % ---------------- required fields ----------------
        useVec = getFieldOr(chanDat,'use',[]);
        if isempty(useVec)
            warning('Missing chanDat.use; skipping file.');
            continue
        end
        useVec = double(useVec(:))==1;
        nBreaths = numel(useVec);
        nUse = sum(useVec);

        behDat = getFieldOr(chanDat,'behDat',[]);
        if isempty(behDat)
            warning('Missing chanDat.behDat; skipping file.');
            continue
        end

        cond = getBehCol(behDat,'condition',nBreaths);
        if sum(isnan(cond))==length(cond)
            cond = getBehColCell(behDat, 'sniffLabel', nBreaths); 
            [~,~,cond] = unique(cond, 'stable'); 
        end
        blen = getBehCol(behDat,'length',nBreaths);
        task = getBehColCell(behDat,'task',nBreaths);
        if sum(cellfun(@(x) strcmp('NA', x), task)) == length(task)
            task = getBehColCell(behDat, 'sniffLabel', nBreaths); 
        end
        nmth = getBehColCell(behDat,'noseMouth',nBreaths);
        warp = getBehColCell(behDat,'warp',nBreaths);

        trial = getFieldOr(chanDat,'trial',[]);
        if isempty(trial) || ~isfield(trial,'rsp')
            warning('Missing chanDat.trial.rsp; skipping file.');
            continue
        end

        tvec = getFieldOr(trial,'tim', getFieldOr(chanDat,'tim',[]));
        rspMat = trial.rsp; % raw (requested)
        [rspMat, tvec] = standardizeTimeByBreath(rspMat, tvec, nBreaths);

        tf = getFieldOr(chanDat,'tf',[]);
        if isempty(tf) || ~isfield(tf,'powZ') || ~isfield(tf,'breathSeg') || ~isfield(tf,'frex')
            warning('Missing chanDat.tf fields; skipping file.');
            continue
        end
        frex = double(tf.frex(:));
        nF = numel(frex);

        powZ = standardize3D(tf.powZ, nBreaths, 50, nF);        % breaths x 50 x frex
        breathSeg = standardize2D(tf.breathSeg, nBreaths, 50);  % breaths x 50

        fooof = getFieldOr(chanDat,'fooof',[]);
        haveFooof = ~isempty(fooof) && isfield(fooof,'aperiodic_fit') && isfield(fooof,'spectra_flat') && isfield(fooof,'gamma_peak_freq');

        reasonElim = getFieldOr(chanDat,'reasonEliminate',[]);
        if ~isempty(reasonElim)
            reasonElim = standardize2D(reasonElim, nBreaths, 4); % breaths x 4
        end

        % ---------------- figure ----------------
        fig = figure('Color','w','Units','pixels','Position',[80 60 1700 1200]);
        if ~visibleFigs, set(fig,'Visible','off'); end

        tl = tiledlayout(fig,4,2,'TileSpacing','compact','Padding','compact');

        titleStr = sprintf('TF/QC summary: %s | sess=%d | task=%s | %s | breaths=%d (QC pass=%d)', ...
            subID, sessNum, taskLab, chanLabel, nBreaths, nUse);
        sgtitle(tl, titleStr, 'FontWeight','bold', 'Interpreter','none');

        % ---------- Panel 1: QC reasons proportion ----------
        ax = nexttile(tl,1); styleAx(ax);
        if isempty(reasonElim)
            text(0.5,0.5,'reasonEliminate missing','HorizontalAlignment','center'); axis off
        else
            counts = sum(reasonElim==1,1,'omitnan');
            props  = counts ./ nBreaths;

            bar(ax, props, 'FaceAlpha',0.9);
            ax.YLim = [0 1];
            ax.XTick = 1:4;
            ax.XTickLabel = {'bad breath','template fail','bad EEG','blink'};
            ax.XTickLabelRotation = 20;
            ylabel(ax,'Proportion');
            title(ax,'QC exclusions (by reason)');
            grid(ax,'on'); ax.GridAlpha = 0.15;

            for k = 1:4
                text(ax, k, props(k)+0.03, sprintf('%d', counts(k)), ...
                    'HorizontalAlignment','center','FontSize',9);
            end
            elimTotal = sum(~useVec);
            text(ax, 0.02, 0.02, sprintf('Total eliminated = %d (overlaps possible)', elimTotal), ...
                'Units','normalized','FontSize',9,'Color',[0.25 0.25 0.25]);
        end

        % ---------- Panel 2: QC-passed breaths per condition ----------
        ax = nexttile(tl,2); styleAx(ax);
        uCond = unique(cond(isfinite(cond)));
        uCond = uCond(:)';
        nCond = numel(uCond);

        countsPass = zeros(1,nCond);
        labelsCond = strings(1,nCond);

        for k = 1:nCond
            c = uCond(k);
            idx = (cond==c) & isfinite(cond);
            countsPass(k) = sum(useVec & idx);
            j = find(idx,1,'first');

            if isnan(warp{j})
                curWarp = ''; 
            else
                curWarp = warp{j}; 
            end
            if strcmp(task{j}, 'shadow')
                if strcmp(chanDat.behDat.shadowFile{j}, 'focusedResp')
                    curTask = 'FocShadow';
                else
                    curTask = 'AudShadow';
                end
            else
                curTask = task{j}; 
            end
            if ~isempty(j)
                labelsCond(k) = [num2str(c) ...
                    curTask nmth{j}, ...
                   curWarp];
            else
                labelsCond(k) = sprintf('c%d', round(c));
            end
        end

        bar(ax, countsPass, 'FaceAlpha',0.9);
        ax.XTick = 1:nCond;
        ax.XTickLabel = labelsCond;
        ax.XTickLabelRotation = 25;
        ylabel(ax,'# breaths');
        title(ax,'QC-passed breaths by condition');
        grid(ax,'on'); ax.GridAlpha = 0.15;

        for k = 1:nCond
            text(ax, k, countsPass(k)+max(1,0.03*max(countsPass)), sprintf('%d', countsPass(k)), ...
                'HorizontalAlignment','center','FontSize',9);
        end
        if nCond > 8
            ax.FontSize = 8; % keep it readable
        end

        % ---------- Panel 3: Breath length density by condition ----------
        ax = nexttile(tl,3); styleAx(ax);
        hold(ax,'on');
        title(ax,'Breath length density (QC-passed)');
        xlabel(ax,'Breath length (s)'); ylabel(ax,'Density');
        grid(ax,'on'); ax.GridAlpha = 0.15;

        legendEntries = {};
        plottedAny = false;
        for k = 1:nCond
            c = uCond(k);
            idx = useVec & (cond==c) & isfinite(blen) & blen>0;
            if sum(idx) < minBreathsPerCond, continue; end
            x = double(blen(idx));

            plotDensity(ax, x);
            legendEntries{end+1} = sprintf('c%d', round(c)); 
            plottedAny = true;
        end

        if plottedAny
            legend(ax, legendEntries, 'Location','best', 'Box','off', 'Interpreter','none');
        else
            cla(ax);
            text(0.5,0.5,'No conditions with enough QC breaths','HorizontalAlignment','center');
            axis(ax,'off');
        end
        hold(ax,'off');

        % ---------- Panel 4: Mean respiration trace by condition ----------
        ax = nexttile(tl,4); styleAx(ax);
        hold(ax,'on');
        title(ax,'Mean respiration trace by condition (QC-passed)');
        xlabel(ax,'Time (s)'); ylabel(ax,'Resp (a.u.)');
        grid(ax,'on'); ax.GridAlpha = 0.15;

        legendEntries = {};
        plottedAny = false;
        for k = 1:nCond
            c = uCond(k);
            idx = useVec & (cond==c);
            if sum(idx) < minBreathsPerCond, continue; end
            m = mean(rspMat(:,idx), 2, 'omitnan');
            plot(ax, tvec, m, 'LineWidth',1.6);
            legendEntries{end+1} = sprintf('c%d', round(c)); 
            plottedAny = true;
        end

        if plottedAny
            legend(ax, legendEntries, 'Location','northeast', 'Box','off', 'Interpreter','none');
        else
            cla(ax);
            text(0.5,0.5,'No conditions with enough QC breaths','HorizontalAlignment','center');
            axis(ax,'off');
        end
        hold(ax,'off');

        % ---------- Panel 5: Mean TF heatmap (powZ) + mean breathSeg ----------
        ax1 = nexttile(tl,5); styleAx(ax1);
        title(ax1,'Mean TF (powZ), QC-passed (all conditions)');
        xlabel(ax1,'Phase-sampled index (1..50)');
        ylabel(ax1,'Frequency (Hz, log)');

        meanPowZ = squeeze(mean(powZ(useVec,:,:), 1, 'omitnan')); % [50 x nF]
        meanPowZ = meanPowZ'; % [nF x 50]
        x = 1:50;

        ylog = log10(frex);
        imagesc(ax1, x, [], meanPowZ);
        axis(ax1,'xy');
        colormap(ax1, parula);
        cb = colorbar(ax1);
        cb.Label.String = 'powZ';
        if ~isempty(tfCLim), caxis(ax1, tfCLim); end

        % log-frequency ticks
        yticks([25:50:300])
        yticklabels(round(frex(25:50:300)))
        ylabel('Frequency (Hz)')

        % epoch boundaries
        hold(ax1,'on');
        for e = phaseEdges
            xline(ax1, e, '-', 'LineWidth',1.0);
        end

        % overlay mean phase-sampled respiration (breathSeg)
        meanSeg = mean(breathSeg(useVec,:), 1, 'omitnan'); % 1 x 50
        yyaxis right
        plot(x, meanSeg, 'k-', 'LineWidth',1.8);
       
        ylabel('Mean resp (a.u.)');
      
        hold(ax1,'off');

        % ---------- Panel 6: Mean aperiodic_fit spectra per epoch ----------
        ax = nexttile(tl,6); styleAx(ax);
        title(ax,'Mean aperiodic\_fit spectra (QC-passed; all conditions)');
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power (a.u.)');
        grid(ax,'on'); ax.GridAlpha = 0.15;
        hold(ax,'on');

        if haveFooof
            aper = standardize3D(fooof.aperiodic_fit, nBreaths, 5, nF); % breaths x 5 x frex
            mAper = squeeze(mean(aper(useVec,:,:), 1, 'omitnan'));      % 5 x nF
            mAper(mAper<=0) = NaN;

            for e = 1:5
                plotSpectrum(ax, frex, mAper(e,:), spectraScale, 1.6);
            end
            legend(ax, epochNames, 'Location','best', 'Box','off', 'Interpreter','none');
           
        else
            cla(ax); text(0.5,0.5,'fooof.aperiodic_fit missing','HorizontalAlignment','center'); axis(ax,'off');
        end
        hold(ax,'off');

        % ---------- Panel 7: Mean spectra_flat spectra per epoch ----------
        ax = nexttile(tl,7); styleAx(ax);
        title(ax,'Mean spectra\_flat spectra (QC-passed; all conditions)');
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Power (a.u.)');
        grid(ax,'on'); ax.GridAlpha = 0.15;
        hold(ax,'on');

        if haveFooof
            flat = standardize3D(fooof.spectra_flat, nBreaths, 5, nF);
            mFlat = squeeze(mean(flat(useVec,:,:), 1, 'omitnan'));
            mFlat(mFlat<=0) = NaN;
            frexMask = frex >=2; 
            for e = 1:5
                plotSpectrum(ax, frex(frexMask), mFlat(e,frexMask), spectraScale, 1.6);
            end
            legend(ax, epochNames, 'Location','best', 'Box','off', 'Interpreter','none');
           
        else
            cla(ax); text(0.5,0.5,'fooof.spectra_flat missing','HorizontalAlignment','center'); axis(ax,'off');
        end
        hold(ax,'off');

        % ---------- Panel 8: Gamma peak frequency densities by epoch ----------
        ax = nexttile(tl,8); styleAx(ax);
        title(ax,'Gamma peak frequency density (QC-passed)');
        xlabel(ax,'Peak frequency (Hz)'); ylabel(ax,'Density');
        grid(ax,'on'); ax.GridAlpha = 0.15;
        hold(ax,'on');

        if haveFooof
            gpf = standardize2D(fooof.gamma_peak_freq, nBreaths, 5); % breaths x 5

            legendEntries = {};
            plottedAny = false;
            for e = 1:5
                xg = double(gpf(useVec,e));
                xg = xg(isfinite(xg));
                if numel(xg) < minBreathsPerCond, continue; end

                plotDensity(ax, xg);
                legendEntries{end+1} = sprintf('%s (n=%d)', epochNames{e}, numel(xg)); %#ok<AGROW>
                plottedAny = true;
            end

            if plottedAny
                legend(ax, legendEntries, 'Location','best', 'Box','off', 'Interpreter','none');
                xlim(ax,[25 60]); % consistent with your gamma peak search band
            else
                cla(ax); text(0.5,0.5,'No gamma peaks with enough QC breaths','HorizontalAlignment','center'); axis(ax,'off');
            end
        else
            cla(ax); text(0.5,0.5,'fooof.gamma_peak_freq missing','HorizontalAlignment','center'); axis(ax,'off');
        end
        hold(ax,'off');

        % ---------------- export ----------------
        try
            exportgraphics(fig, outFile, 'Resolution', dpi);
        catch
            print(fig, outFile, '-djpeg', sprintf('-r%d', dpi));
        end
        close(fig);

        fprintf('  -> saved %s\n', outName);
end %skip already done files! 
    catch ME
        fprintf('  !! ERROR on %s\n', files(ii).name);
        fprintf('     %s\n', ME.message);
        try, close(gcf); end %#ok<TRYNC>
        continue
    end
end

fprintf('Done.\n');

end

%% ================= HELPER FUNCTIONS =================

function v = getFieldOr(S, fld, defaultVal)
if isstruct(S) && isfield(S,fld)
    v = S.(fld);
else
    v = defaultVal;
end
end

function s = safeStr(x)
try
    if isstring(x), s = char(x); return; end
    if ischar(x), s = x; return; end
    if iscell(x) && ~isempty(x), s = safeStr(x{1}); return; end
    if isnumeric(x) && isscalar(x), s = num2str(x); return; end
    s = char(string(x));
catch
    s = 'NA';
end
if isempty(s), s = 'NA'; end
end

function lab = getChanLabel(chanDat, chi)
lab = '';
try
    labels = getFieldOr(chanDat,'labels',{});
    if isempty(labels)
        lab = sprintf('chi%03d', round(chi));
        return
    end
    if isfinite(chi) && chi>=1 && chi<=numel(labels)
        lab = safeStr(labels{round(chi)});
    else
        lab = sprintf('chi%03d', round(chi));
    end
catch
    if isfinite(chi)
        lab = sprintf('chi%03d', round(chi));
    else
        lab = 'chan';
    end
end
if isempty(lab), lab = 'chan'; end
end

function s = sanitizeForFilename(s)
% Keep letters, numbers, underscore, dash; convert others to dash; trim repeats.
s = safeStr(s);
s = strtrim(s);
s = regexprep(s, '\s+', '_');
s = regexprep(s, '[^A-Za-z0-9_\-]', '-');
s = regexprep(s, '-{2,}', '-');
s = regexprep(s, '_{2,}', '_');
if isempty(s), s = 'chan'; end
end

function col = getBehCol(behDat, varName, nBreaths)
col = nan(nBreaths,1);
try
    if istable(behDat) && any(strcmp(behDat.Properties.VariableNames,varName))
        tmp = behDat.(varName);
        col = double(tmp(:));
    elseif isstruct(behDat) && isfield(behDat,varName)
        tmp = behDat.(varName);
        col = double(tmp(:));
    end
catch
end
col = padOrTruncate(col, nBreaths, NaN);
end

function col = getBehColCell(behDat, varName, nBreaths)
col = cell(nBreaths,1); [col{:}] = deal('NA');
try
    if istable(behDat) && any(strcmp(behDat.Properties.VariableNames,varName))
        tmp = behDat.(varName);
    elseif isstruct(behDat) && isfield(behDat,varName)
        tmp = behDat.(varName);
    else
        tmp = [];
    end
    if isempty(tmp), return; end

    tmp = tmp(:);
    for i = 1:min(nBreaths,numel(tmp))
        if iscell(tmp)
            col{i} = safeStr(tmp{i});
        else
            col{i} = safeStr(tmp(i));
        end
    end
catch
end
end

function A = padOrTruncate(A, n, fillVal)
A = A(:);
if numel(A) < n
    A = [A; fillVal*ones(n-numel(A),1)];
elseif numel(A) > n
    A = A(1:n);
end
end

function [X, tvec] = standardizeTimeByBreath(X, tvec, nBreaths)
if isempty(tvec), tvec = (0:size(X,1)-1)'; end
tvec = double(tvec(:));
sz = size(X);

% Expect X is time x breaths or breaths x time; infer by nBreaths
if sz(2) == nBreaths
    X = double(X);
elseif sz(1) == nBreaths
    X = double(X');
else
    % fallback: match by tvec length
    if sz(1) == numel(tvec)
        X = double(X);
    elseif sz(2) == numel(tvec)
        X = double(X');
    else
        error('Cannot infer trial matrix orientation: size=%s, nBreaths=%d', mat2str(sz), nBreaths);
    end
end

% fix tvec length if mismatch
if size(X,1) ~= numel(tvec)
    tvec = (0:size(X,1)-1)';
end
end

function A = standardize2D(A, nRow, nCol)
A = double(A);
sz = size(A);
if isequal(sz, [nRow nCol])
    return
elseif isequal(sz, [nCol nRow])
    A = A.';
    return
else
    if numel(A) == nRow*nCol
        A = reshape(A, [nRow nCol]);
        return
    end
    error('Cannot standardize2D: got %s, need [%d %d]', mat2str(sz), nRow, nCol);
end
end

function A = standardize3D(A, d1, d2, d3)
A = double(A);
sz = size(A);
if numel(sz) ~= 3
    error('standardize3D expects 3D array, got size=%s', mat2str(sz));
end
permsIdx = perms(1:3);
for i = 1:size(permsIdx,1)
    p = permsIdx(i,:);
    if sz(p(1))==d1 && sz(p(2))==d2 && sz(p(3))==d3
        A = permute(A, p);
        return
    end
end
error('Cannot standardize3D: got %s, need permutation of [%d %d %d]', mat2str(sz), d1, d2, d3);
end

function styleAx(ax)
set(ax, 'FontName','Arial', 'FontSize',10, 'LineWidth',1);
box(ax,'off');
end

function plotSpectrum(ax, frex, y, mode, lw)
y = double(y(:))';
frex = double(frex(:))';

switch lower(mode)
    case 'xlinear_ylog'
        semilogy(ax, frex, y, 'LineWidth', lw);
        ax.XScale = 'linear';
        ax.YScale = 'log';
    case 'linear'
        plot(ax, frex, y, 'LineWidth', lw);
        ax.XScale = 'linear';
        ax.YScale = 'linear';
    otherwise
        % fallback
        semilogy(ax, frex, y, 'LineWidth', lw);
        ax.XScale = 'linear';
        ax.YScale = 'log';
end
end


function plotDensity(ax, x)
x = double(x(:));
x = x(isfinite(x));
if isempty(x), return; end
if exist('ksdensity','file') == 2
    [f,xi] = ksdensity(x);
    plot(ax, xi, f, 'LineWidth',1.5);
else
    histogram(ax, x, 'Normalization','pdf', 'DisplayStyle','stairs', 'LineWidth',1.5);
end
end
