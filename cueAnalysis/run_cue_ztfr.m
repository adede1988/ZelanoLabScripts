function run_cue_ztfr(sessFilter)
% RUN_CUE_ZTFR  (analysis3) Per session: detect sharp-deflection noise trials,
%   draw the red-flagged single-trial raw plot, then compute the bestMac
%   bootstrap-z spectrograms (trialStart + finalOnset, shared trialStart
%   baseline, respiration overlay) on the NOISE-FREE trials. Overwrites the
%   per-subject figures/.mat and writes the noise counts to cueTask_fooof_summary.csv.

    if nargin < 1, sessFilter = []; end
    cue_init_paths(); L = labPaths();
    groupDir = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs';
    csvPath  = fullfile(groupDir, 'cueTask_fooof_summary.csv');
    dispWin  = [-1000 3000]; sep = 50;

    T = cue_session_table(false); T = T(T.onDisk, :);
    if ~isempty(sessFilter), T = T(ismember(T.sessID, string(sessFilter)), :); end

    % augment FOOOF CSV with the analysis3 columns; null the obsolete z-QC ones
    F = readtable(csvPath); F.sessID = string(F.sessID);
    newCols = {'nTrials','nNoiseTrials','nOOB_finalOnset','nFinal_trialStart','nFinal_finalOnset'};
    for c = 1:numel(newCols)
        if ~ismember(newCols{c}, F.Properties.VariableNames), F.(newCols{c}) = nan(height(F),1); end
    end
    for old = {'nDropBaseZ10','nDropPostZ15_trialStart','nDropPostZ15_finalOnset','nOOB_trialStart','nDrop_noBaseline_finalOnset'}
        if ismember(old{1}, F.Properties.VariableNames), F.(old{1})(:) = NaN; end
    end

    for i = 1:height(T)
        id = char(T.sessID(i)); fp = char(T.path(i));
        fprintf('\n== %d/%d %s ==\n', i, height(T), id);
        try
            Sv = load(fp); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;
            if ~isfield(od,'bestMac') || isempty(od.bestMac), fprintf('  no bestMac -> skip\n'); clear od; continue; end
            labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
            ci = find(strcmp(od.bestMac, labs), 1); if isempty(ci), clear od; continue; end
            sig = double(od.data(ci, :)); fs = od.fs;

            isR = cellfun(@(x) contains(x,'rsp'), labs);
            rspAll = od.data(isR, :);
            ridx = od.rspIDX; if isempty(ridx) || ridx < 1 || ridx > size(rspAll,1), ridx = 1; end
            rsp = double(rspAll(ridx, :)) .* od.rspFlip;

            figDir = fullfile(L.figPath, id, 'cueTask'); if ~isfolder(figDir), mkdir(figDir); end

            % ---- noise detection on ALL trialStart trials + single-trial plot ----
            ts = od.TTL.trialStart;
            NT = cue_noise_trials(sig, fs, ts);
            cue_plot_singletrial(NT, id, od.bestMac, sep, fullfile(figDir, 'singleTrialRawMac.png'));
            nNoise = sum(NT.noisy(NT.ok));

            % record noise counts NOW so they survive even if spectrograms are skipped
            sel = (F.sessID == string(id)) & (F.isBestMac == 1);
            F.nTrials(sel) = numel(ts); F.nNoiseTrials(sel) = nNoise;
            F.nFinal_trialStart(sel) = 0; F.nFinal_finalOnset(sel) = 0; F.nOOB_finalOnset(sel) = NaN;
            writetable(F, csvPath);

            matFile = fullfile(figDir, [id '_cue_bestMac_TFR.mat']);

            % ---- paired trialStart/finalOnset (behDat), drop NOISY + invalid ----
            bd = od.behDat;
            tsV = nan(height(bd),1); foV = nan(height(bd),1); noiseFlag = false(height(bd),1);
            for j = 1:height(bd)
                k = bd.n(j);
                if k >= 1 && k <= numel(ts)
                    tsV(j) = ts(k); foV(j) = bd.finalOnset(j);
                    if k <= numel(NT.noisy), noiseFlag(j) = NT.noisy(k); end
                end
            end
            keep = isfinite(tsV) & isfinite(foV) & ~noiseFlag;
            tsV = tsV(keep); foV = foV(keep);

            R = cue_ztfr_pair(sig, rsp, fs, tsV, foV);
            if isempty(R.trialStart)
                % session too noisy: drop the stale TFR .mat so the group step excludes it
                if isfile(matFile), delete(matFile); end
                fprintf('  too few clean trials -> spectrograms skipped, EXCLUDED (%d/%d noise)\n', nNoise, sum(NT.ok));
                clear od; continue;
            end

            a = prctile(abs(R.trialStart.map(:)), 98); if ~isfinite(a) || a <= 0, a = 20; end
            clim = [-a a];
            cue_plot_ztfr(R.trialStart.map, R.trialStart.times, R.freqs, clim, ...
                R.trialStart.resp, R.trialStart.respT, dispWin, ...
                sprintf('%s  %s  trialStart-locked (z, n=%d)', id, od.bestMac, R.trialStart.nFinal), ...
                fullfile(figDir, [id '_' od.bestMac '_TFR_trialStart.png']));
            if ~isempty(R.finalOnset)
                cue_plot_ztfr(R.finalOnset.map, R.finalOnset.times, R.freqs, clim, ...
                    R.finalOnset.resp, R.finalOnset.respT, dispWin, ...
                    sprintf('%s  %s  finalOnset-locked (z, n=%d)', id, od.bestMac, R.finalOnset.nFinal), ...
                    fullfile(figDir, [id '_' od.bestMac '_TFR_finalOnset.png']));
            end

            tfrOut = struct('sessID', id, 'subID', char(T.subID(i)), 'group', char(T.group(i)), ...
                'type', char(T.type(i)), 'sessNum', T.sessNum(i), 'bestMac', od.bestMac, ...
                'dispWin', dispWin, 'clim', clim, 'freqs', R.freqs, ...
                'trialStart', R.trialStart, 'finalOnset', R.finalOnset, 'qc', R.qc, ...
                'nNoiseTrials', nNoise); %#ok<NASGU>
            save(fullfile(figDir, [id '_cue_bestMac_TFR.mat']), 'tfrOut', '-v7');

            F.nOOB_finalOnset(sel) = R.qc.nOOB_finalOnset;
            F.nFinal_trialStart(sel) = R.qc.nFinal_trialStart;
            F.nFinal_finalOnset(sel) = R.qc.nFinal_finalOnset;
            writetable(F, csvPath);

            fprintf('  noise=%d/%d  nFinal ts=%d fo=%d  clim=%.1f\n', ...
                nNoise, sum(NT.ok), R.qc.nFinal_trialStart, R.qc.nFinal_finalOnset, a);
            clear od sig rsp R tfrOut NT;
        catch ME
            fprintf('  FAILED: %s\n', ME.message); try, clear od; catch, end
        end
    end
    fprintf('\nanalysis3 z-TFR pass done. CSV: %s\n', csvPath);
end
