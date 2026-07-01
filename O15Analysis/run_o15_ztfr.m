function run_o15_ztfr(sessFilter)
% RUN_O15_ZTFR  Per O15 session: split sniffs into start/free/confirm, flag
%   RELATIVE sharp-deflection noise (cue_noise_trials, robust-z zd>K with K from
%   o15_noise_K), draw the red-flagged single-trial raw plot for the START
%   sniffs, then compute the bestMac
%   baseline-z spectrograms for ALL THREE sniff types (o15_ztfr_multi) against
%   ONE shared baseline taken before the first start sniff. Saves the per-subject
%   PNGs + numeric .mat and writes the noise/QC columns to O15Task_fooof_summary.csv.
%
%   Lockings: start / free / confirm (behDat.finalOnset by sniffLabel).
%   Reuses cue_noise_trials, cue_plot_singletrial, cue_plot_ztfr unchanged.

    if nargin < 1, sessFilter = []; end
    o15_init_paths(); L = labPaths();
    groupDir = getenv('O15_GROUPDIR');
    if isempty(groupDir), groupDir = fullfile(labPaths().figPath, 'groupStatFigs'); end
    csvPath  = fullfile(groupDir, 'O15Task_fooof_summary.csv');
    dispWin  = [-1000 3000]; sep = 50;
    types = {'start','free','confirm'};

    T = o15_session_table(false); T = T(T.onDisk & T.fresh, :);
    if ~isempty(sessFilter), T = T(ismember(T.sessID, string(sessFilter)), :); end

    % augment FOOOF CSV with the O15 QC columns
    F = readtable(csvPath); F.sessID = string(F.sessID);
    qcCols = {'nStart','nFree','nConfirm','nNoiseStart','nNoiseFree','nNoiseConfirm', ...
              'nFinal_start','nFinal_free','nFinal_confirm','baseFrames','baseDurS'};
    for c = 1:numel(qcCols)
        if ~ismember(qcCols{c}, F.Properties.VariableNames), F.(qcCols{c}) = nan(height(F),1); end
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

            figDir = fullfile(L.figPath, id, 'O15'); if ~isfolder(figDir), mkdir(figDir); end
            sel = (F.sessID == string(id)) & (F.isBestMac == 1);

            % ---- split sniffs by label + one noise pass over all sniff onsets ----
            bd = od.behDat;
            sl = string(bd.sniffLabel); fo = bd.finalOnset;
            Knoise = o15_noise_K();                  % relative threshold (robust-SDs), single source
            NTall = cue_noise_trials(sig, fs, fo, [], Knoise);   % per-sniff RELATIVE noise flag (row-aligned)

            onsets = struct(); nNoise = struct();
            for ti = 1:numel(types)
                tp = types{ti};
                rows = find(sl == tp);
                good = rows(NTall.ok(rows) & ~NTall.noisy(rows));
                onsets.(tp) = fo(good);
                nNoise.(tp) = sum(NTall.ok(rows) & NTall.noisy(rows));
                F.(['n' upper(tp(1)) tp(2:end)])(sel) = numel(rows);
                F.(['nNoise' upper(tp(1)) tp(2:end)])(sel) = nNoise.(tp);
            end
            writetable(F, csvPath);   % counts survive even if spectrograms skipped

            % ---- single-trial raw QC plot for START sniffs (red = noise) ----
            NTstart = cue_noise_trials(sig, fs, fo(sl=="start"), [], Knoise);
            cue_plot_singletrial(NTstart, id, od.bestMac, sep, fullfile(figDir, 'singleTrialRawMac_start.png'));

            matFile = fullfile(figDir, [id '_O15_bestMac_TFR.mat']);

            % ---- baseline-z spectrograms for the three sniff types ----
            R = o15_ztfr_multi(sig, rsp, fs, onsets);
            if isempty(R.freqs) || isempty(R.baseline) || R.baseline.nFrames == 0
                if isfile(matFile), delete(matFile); end
                fprintf('  no usable baseline / too few clean sniffs -> EXCLUDED\n'); clear od; continue;
            end

            % shared color scale across the three maps (within-session comparability)
            allmaps = [];
            for ti = 1:numel(types), if ~isempty(R.(types{ti})), allmaps = [allmaps; abs(R.(types{ti}).map(:))]; end, end %#ok<AGROW>
            a = prctile(allmaps, 98); if ~isfinite(a) || a <= 0, a = 5; end
            clim = [-a a];

            for ti = 1:numel(types)
                tp = types{ti}; S = R.(tp);
                if isempty(S), continue; end
                cue_plot_ztfr(S.map, S.times, S.freqs, clim, S.resp, S.respT, dispWin, ...
                    sprintf('%s  %s  %s-sniff (z, n=%d)', id, od.bestMac, tp, S.nFinal), ...
                    fullfile(figDir, [id '_' od.bestMac '_O15TFR_' tp '.png']));
                F.(['nFinal_' tp])(sel) = S.nFinal;
            end
            F.baseFrames(sel) = R.baseline.nFrames;
            F.baseDurS(sel)   = R.baseline.durS;
            writetable(F, csvPath);

            tfrOut = struct('sessID', id, 'subID', char(T.subID(i)), 'group', char(T.group(i)), ...
                'type', char(T.type(i)), 'sessNum', T.sessNum(i), 'bestMac', od.bestMac, ...
                'dispWin', dispWin, 'clim', clim, 'freqs', R.freqs, 'baseline', R.baseline, ...
                'start', R.start, 'free', R.free, 'confirm', R.confirm, 'qc', R.qc); %#ok<NASGU>
            save(matFile, 'tfrOut', '-v7');

            fprintf('  base=%d frames (%.1f s)  nFinal start/free/confirm = %d/%d/%d  clim=%.2f\n', ...
                R.baseline.nFrames, R.baseline.durS, gv(R.start), gv(R.free), gv(R.confirm), a);
            clear od sig rsp R tfrOut NTall NTstart;
        catch ME
            fprintf('  FAILED: %s\n', ME.message);
            fprintf('  @ %s line %d\n', ME.stack(1).name, ME.stack(1).line);
            try, clear od; catch, end
        end
    end
    fprintf('\nO15 baseline-z TFR pass done. CSV: %s\n', csvPath);
end

function n = gv(S), if isempty(S), n = 0; else, n = S.nFinal; end, end
