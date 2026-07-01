function run_thresh_ztfr(sessFilter)
% RUN_THRESH_ZTFR  Per session: detect sharp-deflection noise trials on the
%   TTL.start-locked bestMac epochs, draw the red-flagged single-trial raw plot,
%   then compute the bestMac bootstrap-z spectrograms split by ODOR CONDITION
%   (low / med / air), all finalOnset-locked, each z-scored against its own trials'
%   pre-`start` baseline (thresh_ztfr_conds), with the mean respiration overlaid.
%   Overwrites the per-subject figures/.mat and writes noise/QC counts to
%   threshTask_fooof_summary.csv.
%
%   The single scientific change from run_cue_ztfr: the finalOnset response is
%   split three ways by behDat.type instead of paired against a trialStart locking.

    if nargin < 1, sessFilter = []; end
    thresh_init_paths(); L = labPaths();
    groupDir = getenv('THRESH_GROUPDIR');                    % per-job override (parallel runs)
    if isempty(groupDir), groupDir = fullfile(L.figPath, 'groupStatFigs'); end
    csvPath  = fullfile(groupDir, 'threshTask_fooof_summary.csv');
    dispWin  = [-1000 3000]; sep = 50;
    conds = {'low','med','air'};

    T = thresh_session_table(false); T = T(T.onDisk, :);
    if ~isempty(sessFilter), T = T(ismember(T.sessID, string(sessFilter)), :); end

    % augment FOOOF CSV with the condition-split QC columns
    F = readtable(csvPath); F.sessID = string(F.sessID);
    newCols = {'nTrials','nNoiseTrials','nLow','nMed','nAir', ...
               'nFinal_low','nFinal_med','nFinal_air'};
    for c = 1:numel(newCols)
        if ~ismember(newCols{c}, F.Properties.VariableNames), F.(newCols{c}) = nan(height(F),1); end
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

            figDir = fullfile(L.figPath, id, 'threshTask'); if ~isfolder(figDir), mkdir(figDir); end

            % ---- noise detection on ALL TTL.start trials + single-trial plot ----
            ts = od.TTL.start;
            NT = thresh_noise_trials(sig, fs, ts);
            thresh_plot_singletrial(NT, id, od.bestMac, sep, fullfile(figDir, 'singleTrialRawMac.png'));
            nNoise = sum(NT.noisy(NT.ok));

            % record noise counts NOW so they survive even if spectrograms are skipped
            sel = (F.sessID == string(id)) & (F.isBestMac == 1);
            F.nTrials(sel) = numel(ts); F.nNoiseTrials(sel) = nNoise;
            F.nLow(sel) = 0; F.nMed(sel) = 0; F.nAir(sel) = 0;
            F.nFinal_low(sel) = 0; F.nFinal_med(sel) = 0; F.nFinal_air(sel) = 0;
            writetable(F, csvPath);

            matFile = fullfile(figDir, [id '_thresh_bestMac_TFR.mat']);

            % ---- build paired (start, finalOnset) vectors per condition from behDat,
            %      dropping NOISY (by trial anchor k) + invalid ----
            bd = od.behDat;
            startByCond = struct('low',[],'med',[],'air',[]);
            foByCond    = struct('low',[],'med',[],'air',[]);
            typ = string(bd.type);
            for j = 1:height(bd)
                k = bd.n(j);
                if ~(k >= 1 && k <= numel(ts)), continue; end
                s = ts(k); f = bd.finalOnset(j); cond = char(typ(j));
                if ~ismember(cond, conds), continue; end
                if ~(isfinite(s) && isfinite(f)), continue; end
                if k <= numel(NT.noisy) && NT.noisy(k), continue; end   % drop noisy trial
                startByCond.(cond)(end+1,1) = s;
                foByCond.(cond)(end+1,1)    = f;
            end
            nLow = numel(foByCond.low); nMed = numel(foByCond.med); nAir = numel(foByCond.air);
            F.nLow(sel) = nLow; F.nMed(sel) = nMed; F.nAir(sel) = nAir;
            writetable(F, csvPath);

            % ---- per-condition bootstrap-z spectrograms ----
            R = thresh_ztfr_conds(sig, rsp, fs, startByCond, foByCond);
            if isempty(R.freqs)
                % session too noisy / too few trials in every condition: drop stale mat
                if isfile(matFile), delete(matFile); end
                fprintf('  too few clean trials in all conditions -> spectrograms skipped, EXCLUDED (%d/%d noise)\n', nNoise, sum(NT.ok));
                clear od; continue;
            end

            % shared color scale across the three condition maps (comparability)
            allmaps = [];
            for cc = 1:numel(conds)
                S = R.(conds{cc});
                if ~isempty(S), allmaps = [allmaps; abs(S.map(:))]; end %#ok<AGROW>
            end
            a = prctile(allmaps, 98); if ~isfinite(a) || a <= 0, a = 20; end
            clim = [-a a];

            for cc = 1:numel(conds)
                cond = conds{cc}; S = R.(cond);
                if isempty(S), continue; end
                thresh_plot_ztfr(S.map, S.times, S.freqs, clim, S.resp, S.respT, dispWin, ...
                    sprintf('%s  %s  %s-odor  finalOnset-locked (z, n=%d)', id, od.bestMac, cond, S.nFinal), ...
                    fullfile(figDir, [id '_' od.bestMac '_TFR_' cond '.png']));
                F.(['nFinal_' cond])(sel) = S.nFinal;
            end
            writetable(F, csvPath);

            tfrOut = struct('sessID', id, 'subID', char(T.subID(i)), 'group', char(T.group(i)), ...
                'type', char(T.type(i)), 'sessNum', T.sessNum(i), 'bestMac', od.bestMac, ...
                'dispWin', dispWin, 'clim', clim, 'freqs', R.freqs, ...
                'low', R.low, 'med', R.med, 'air', R.air, 'qc', R.qc, ...
                'nNoiseTrials', nNoise); %#ok<NASGU>
            save(matFile, 'tfrOut', '-v7');

            fprintf('  noise=%d/%d  nFinal low/med/air = %d/%d/%d  clim=%.1f\n', ...
                nNoise, sum(NT.ok), R.qc.low.nFinal, R.qc.med.nFinal, R.qc.air.nFinal, a);
            clear od sig rsp R tfrOut NT;
        catch ME
            fprintf('  FAILED: %s\n', ME.message);
            if ~isempty(ME.stack), fprintf('  @ %s line %d\n', ME.stack(1).name, ME.stack(1).line); end
            try, clear od; catch, end
        end
    end
    fprintf('\ncondition-split z-TFR pass done. CSV: %s\n', csvPath);
end
