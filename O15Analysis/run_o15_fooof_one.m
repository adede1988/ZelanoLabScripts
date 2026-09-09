function run_o15_fooof_one(sessIDs, doSaveBack)
% RUN_O15_FOOOF_ONE  FOOOF + bestMac for specific O15 session(s) and APPEND their
%   rows to the existing O15Task_fooof_summary.csv (does NOT overwrite other
%   sessions). Use to integrate a newly-completed participant before run_o15_ztfr.
%   Mirrors the FOOOF half of run_o15_fooof_all for the named sessions.

    if nargin < 2 || isempty(doSaveBack), doSaveBack = true; end
    o15_init_paths(); L = labPaths();
    groupDir = getenv('O15_GROUPDIR');
    if isempty(groupDir), groupDir = fullfile(labPaths().figPath, 'groupStatFigs'); end
    csvPath  = fullfile(groupDir, 'O15Task_fooof_summary.csv');
    gammaBand = [30 58];

    T = o15_session_table(false); T = T(T.onDisk & T.fresh, :);
    T = T(ismember(T.sessID, string(sessIDs)), :);
    if isempty(T), error('run_o15_fooof_one:none', 'No fresh on-disk sessions matched.'); end

    hdr = {'subID','sessID','sessNum','type','group','channel', ...
           'gammaPeakDetected','peakGammaFreq','peakGammaPower','flattenedGammaMax', ...
           'aperiodicExponent','aperiodicOffset','fooofR2','isBestMac','selectionMethod'};

    if isfile(csvPath)
        F = readtable(csvPath);
    else
        F = cell2table(cell(0,numel(hdr)), 'VariableNames', hdr);
    end

    for i = 1:height(T)
        id = char(T.sessID(i)); fp = char(T.path(i));
        fprintf('\n== FOOOF %s ==\n', id);
        Sv = load(fp); fn = fieldnames(Sv); od = Sv.(fn{1}); clear Sv;
        R = cue_fooof_macBP(od, gammaBand);
        if R.empty, fprintf('  no macBP channels -> cannot integrate\n'); clear od; continue; end

        figDir = fullfile(L.figPath, id, 'O15'); if ~isfolder(figDir), mkdir(figDir); end
        cue_plot_fooof(R, id, fullfile(figDir, [id '_macBP_fooof_periodic.png']));

        od.bestMac = R.bestMac; od.bestMacMethod = R.selectionMethod;
        if doSaveBack
            tmp = [fp '.tmp']; SS = struct(); SS.(fn{1}) = od;
            save(tmp, '-struct', 'SS', '-v7.3'); movefile(tmp, fp, 'f');
            fprintf('  wrote bestMac=%s into %s\n', R.bestMac, id);
        end

        rows = cell(numel(R.labels), numel(hdr));
        for m = 1:numel(R.labels)
            rows(m,:) = {char(T.subID(i)), id, T.sessNum(i), char(T.type(i)), char(T.group(i)), ...
                R.labels{m}, double(R.gammaDetected(m)), R.peakGammaFreq(m), R.peakGammaPower(m), ...
                R.flatGammaMax(m), R.apExponent(m), R.apOffset(m), R.r2(m), double(m==R.bestIdx), ...
                R.selectionMethod};
        end
        nr = cell2table(rows, 'VariableNames', hdr);

        % align to F's columns (extra QC cols -> NaN), drop any stale rows, append
        if width(F) > 0
            extra = setdiff(F.Properties.VariableNames, nr.Properties.VariableNames, 'stable');
            for v = 1:numel(extra), nr.(extra{v}) = nan(height(nr),1); end
            nr = nr(:, F.Properties.VariableNames);
            F(strcmpi(string(F.sessID), id), :) = [];
        end
        F = [F; nr]; %#ok<AGROW>
        fprintf('  appended %d channel rows (bestMac=%s, %s)\n', numel(R.labels), R.bestMac, R.selectionMethod);
        clear od R;
    end

    writetable(F, csvPath);
    fprintf('\nO15 FOOOF integration done -> %s (%d rows)\n', csvPath, height(F));
end
