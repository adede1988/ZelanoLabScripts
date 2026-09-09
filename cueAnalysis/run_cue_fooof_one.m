function run_cue_fooof_one(sessIDs, doSaveBack)
% RUN_CUE_FOOOF_ONE  FOOOF + bestMac for specific session(s) and APPEND their
%   rows to the existing cueTask_fooof_summary.csv (does NOT overwrite other
%   sessions). Use to integrate a newly-available participant before run_cue_ztfr.
%   Mirrors the FOOOF half of run_cue_tasks123 (Task 1) for the named sessions.

    if nargin < 2 || isempty(doSaveBack), doSaveBack = true; end
    cue_init_paths(); L = labPaths();
    groupDir = getenv('CUE_GROUPDIR');
    if isempty(groupDir), groupDir = fullfile(L.figPath, 'groupStatFigs'); end
    csvPath  = fullfile(groupDir, 'cueTask_fooof_summary.csv');
    gammaBand = [30 58];

    T = cue_session_table(false); T = T(T.onDisk, :);
    T = T(ismember(T.sessID, string(sessIDs)), :);
    if isempty(T), error('run_cue_fooof_one:none', 'No on-disk sessions matched.'); end

    hdr = {'subID','sessID','sessNum','type','group','channel', ...
           'gammaPeakDetected','peakGammaFreq','peakGammaPower','flattenedGammaMax', ...
           'aperiodicExponent','aperiodicOffset','fooofR2','isBestMac','selectionMethod','rejRate'};

    F = readtable(csvPath);   % text cols -> cellstr, numeric -> double

    for i = 1:height(T)
        id = char(T.sessID(i)); fp = char(T.path(i));
        fprintf('\n== FOOOF %s ==\n', id);
        Sv = load(fp); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;
        R = cue_fooof_macBP(od, gammaBand);
        if R.empty, fprintf('  no macBP channels -> cannot integrate\n'); clear od; continue; end

        figDir = fullfile(L.figPath, id, 'cueTask'); if ~isfolder(figDir), mkdir(figDir); end
        cue_plot_fooof(R, id, fullfile(figDir, [id '_macBP_fooof_periodic.png']));

        od.bestMac = R.bestMac; od.bestMacMethod = R.selectionMethod;
        if doSaveBack
            tmp = [fp '.tmp']; SS = struct(); SS.(fn{1}) = od;
            save(tmp, '-struct', 'SS', '-v7.3'); movefile(tmp, fp, 'f');
            fprintf('  wrote bestMac=%s into %s\n', R.bestMac, id);
        end

        % build this session's channel rows
        rows = cell(numel(R.labels), numel(hdr));
        for m = 1:numel(R.labels)
            rows(m,:) = {char(T.subID(i)), id, T.sessNum(i), char(T.type(i)), char(T.group(i)), ...
                R.labels{m}, double(R.gammaDetected(m)), R.peakGammaFreq(m), R.peakGammaPower(m), ...
                R.flatGammaMax(m), R.apExponent(m), R.apOffset(m), R.r2(m), double(m==R.bestIdx), ...
                R.selectionMethod, R.rejRate(m)};
        end
        nr = cell2table(rows, 'VariableNames', hdr);

        % align to F's columns (extra QC cols -> NaN), drop any stale rows, append
        extra = setdiff(F.Properties.VariableNames, nr.Properties.VariableNames, 'stable');
        for v = 1:numel(extra), nr.(extra{v}) = nan(height(nr),1); end
        nr = nr(:, F.Properties.VariableNames);
        F(strcmpi(string(F.sessID), id), :) = [];
        F = [F; nr]; %#ok<AGROW>
        fprintf('  appended %d channel rows (bestMac=%s, %s)\n', numel(R.labels), R.bestMac, R.selectionMethod);
        clear od R;
    end

    writetable(F, csvPath);
    fprintf('\nFOOOF integration done -> %s (%d rows)\n', csvPath, height(F));
end
