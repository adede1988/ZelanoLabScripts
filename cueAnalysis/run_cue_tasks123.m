function run_cue_tasks123(sessFilter, doSaveBack, csvPath)
% RUN_CUE_TASKS123  FOOOF / bestMac stage (analysis "Task 1") for all sessions.
%   Per cue session: FOOOF the macBP channels, pick the bestMac channel, write
%   bestMac into the .mat, save the per-subject FOOOF periodic-spectrum figure,
%   and (over)write the per-channel rows of cueTask_fooof_summary.csv.
%
%   (Legacy name -- it once also produced dB spectrograms, now superseded by
%   run_cue_ztfr; only the FOOOF half remains here.)
%
%   run_cue_tasks123(sessFilter, doSaveBack, csvPath)
%     sessFilter : cellstr/string of sessIDs ([] = all on disk)
%     doSaveBack : true -> write outDat.bestMac into the .mat (safe temp+rename)
%     csvPath    : output CSV ([] = groupStatFigs\cueTask_fooof_summary.csv)
%
%   WARNING: this REWRITES the whole CSV from the processed sessions' rows. To
%   add a single new participant WITHOUT disturbing the others, use
%   run_cue_fooof_one (which appends/replaces just that session's rows).

    if nargin < 1, sessFilter = []; end
    if nargin < 2 || isempty(doSaveBack), doSaveBack = false; end

    cue_init_paths();
    L = labPaths();
    groupDir = getenv('CUE_GROUPDIR');                       % per-job override (parallel runs)
    if isempty(groupDir), groupDir = fullfile(L.figPath, 'groupStatFigs'); end
    if ~isfolder(groupDir), mkdir(groupDir); end
    if nargin < 3 || isempty(csvPath)
        csvPath = fullfile(groupDir, 'cueTask_fooof_summary.csv');
    end
    gammaBand = [30 58];

    T = cue_session_table(false);
    T = T(T.onDisk, :);
    if ~isempty(sessFilter)
        T = T(ismember(T.sessID, string(sessFilter)), :);
    end

    hdr = {'subID','sessID','sessNum','type','group','channel', ...
           'gammaPeakDetected','peakGammaFreq','peakGammaPower','flattenedGammaMax', ...
           'aperiodicExponent','aperiodicOffset','fooofR2','isBestMac','selectionMethod'};
    rows = {};

    for i = 1:height(T)
        id = char(T.sessID(i)); fp = char(T.path(i));
        fprintf('\n==== %d/%d  %s (%s) ====\n', i, height(T), id, char(T.group(i)));
        try
            Sv = load(fp); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;

            R = cue_fooof_macBP(od, gammaBand);
            if R.empty
                fprintf('  no macBP channels -> skip session\n'); clear od R; continue;
            end

            figDir = fullfile(L.figPath, id, 'cueTask');
            if ~isfolder(figDir), mkdir(figDir); end
            cue_plot_fooof(R, id, fullfile(figDir, [id '_macBP_fooof_periodic.png']));

            for m = 1:numel(R.labels)
                rows(end+1,:) = {char(T.subID(i)), id, T.sessNum(i), char(T.type(i)), ...
                    char(T.group(i)), R.labels{m}, double(R.gammaDetected(m)), ...
                    R.peakGammaFreq(m), R.peakGammaPower(m), R.flatGammaMax(m), ...
                    R.apExponent(m), R.apOffset(m), R.r2(m), double(m==R.bestIdx), ...
                    R.selectionMethod}; %#ok<AGROW>
            end

            od.bestMac = R.bestMac; od.bestMacMethod = R.selectionMethod;
            if doSaveBack
                try
                    saveBackBestMac(fp, fn{1}, od);
                    fprintf('  wrote bestMac=%s into %s\n', R.bestMac, id);
                catch ME
                    fprintf('  SAVE-BACK FAILED (%s): %s\n', id, ME.message);
                end
            end

            % flush the CSV after every session so partial progress survives
            if ~isempty(rows), writetable(cell2table(rows, 'VariableNames', hdr), csvPath); end
            fprintf('  bestMac=%s (%s)\n', R.bestMac, R.selectionMethod);
            clear od R;
        catch ME
            fprintf('  SESSION FAILED (%s): %s\n', id, ME.message);
            try, clear od R; catch, end
        end
    end

    fprintf('\nFOOOF stage done. CSV=%s (%d channel rows)\n', csvPath, size(rows,1));
end

% ------------------------------------------------------------------
function saveBackBestMac(fp, varName, od)
    tmp = [fp '.tmp'];
    Sv = struct(); Sv.(varName) = od;
    save(tmp, '-struct', 'Sv', '-v7.3');
    movefile(tmp, fp, 'f');
end
