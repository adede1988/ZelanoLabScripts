function run_o15_fooof_all(sessFilter, doSaveBack, csvPath)
% RUN_O15_FOOOF_ALL  FOOOF / bestMac stage for the O15 task (all fresh sessions).
%   Per O15 session: FOOOF the macBP channels (REUSES cue_fooof_macBP), pick the
%   bestMac channel, write bestMac into the .mat, save the per-subject FOOOF
%   periodic-spectrum figure (cue_plot_fooof), and (over)write the per-channel
%   rows of O15Task_fooof_summary.csv. Sessions with no macBP channels are
%   skipped (O15 EEG/macros are sometimes absent).
%
%   run_o15_fooof_all(sessFilter, doSaveBack, csvPath)
%     sessFilter : cellstr/string of sessIDs ([] = all fresh on-disk)
%     doSaveBack : true -> write outDat.bestMac into the .mat (safe temp+rename)
%     csvPath    : output CSV ([] = groupStatFigs\O15Task_fooof_summary.csv)
%
%   WARNING: REWRITES the whole CSV from the processed sessions' rows. To add a
%   single new participant without disturbing the others, use run_o15_fooof_one.

    if nargin < 1, sessFilter = []; end
    if nargin < 2 || isempty(doSaveBack), doSaveBack = true; end

    o15_init_paths();
    L = labPaths();
    groupDir = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs';
    if ~isfolder(groupDir), mkdir(groupDir); end
    if nargin < 3 || isempty(csvPath)
        csvPath = fullfile(groupDir, 'O15Task_fooof_summary.csv');
    end
    gammaBand = [30 58];

    T = o15_session_table(false);
    T = T(T.onDisk & T.fresh, :);
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

            R = cue_fooof_macBP(od, gammaBand);     % REUSED: generic over any outDat w/ macBP
            if R.empty
                fprintf('  no macBP channels -> skip session\n'); clear od R; continue;
            end

            figDir = fullfile(L.figPath, id, 'O15');
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

            if ~isempty(rows), writetable(cell2table(rows, 'VariableNames', hdr), csvPath); end
            fprintf('  bestMac=%s (%s)\n', R.bestMac, R.selectionMethod);
            clear od R;
        catch ME
            fprintf('  SESSION FAILED (%s): %s\n', id, ME.message);
            try, clear od R; catch, end
        end
    end

    fprintf('\nO15 FOOOF stage done. CSV=%s (%d channel rows)\n', csvPath, size(rows,1));
end

% ------------------------------------------------------------------
function saveBackBestMac(fp, varName, od)
    tmp = [fp '.tmp'];
    Sv = struct(); Sv.(varName) = od;
    save(tmp, '-struct', 'Sv', '-v7.3');
    movefile(tmp, fp, 'f');
end
