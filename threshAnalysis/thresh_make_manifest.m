function thresh_make_manifest()
% THRESH_MAKE_MANIFEST  Regenerate threshTask_data_manifest.csv from the current
%   FOOOF+QC summary CSV + the session table, so the data-availability table
%   always reflects the live analysis state:
%     status = 'ok'                | macBP present, spectrograms made
%              'noMacBP'           | session has no macBP channels (skipped)
%              'excluded(allNoise)'| macBP present but no clean trials in any cond
%   Derived (not hand-written) so it cannot drift from the rest of the report.

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo, 'threshAnalysis'));
    L = labPaths();
    groupDir = getenv('THRESH_GROUPDIR');
    if isempty(groupDir), groupDir = fullfile(L.figPath, 'groupStatFigs'); end
    csvPath  = fullfile(groupDir, 'threshTask_fooof_summary.csv');
    manPath  = fullfile(groupDir, 'threshTask_data_manifest.csv');

    T = thresh_session_table(false); T = T(T.onDisk, :);
    F = readtable(csvPath); Fsess = string(F.sessID);

    getcol = @(sub, name, bi) local_getcol(sub, name, bi);

    n = height(T);
    rows = cell(n, 12);
    for i = 1:n
        id = char(T.sessID(i));
        sel = strcmpi(Fsess, id);
        if ~any(sel)
            status = 'noMacBP'; nMac = 0; bestMac = ''; nNoise = NaN; nFl = NaN; nFm = NaN; nFa = NaN;
        else
            sub = F(sel, :);
            nMac = height(sub);
            bi = find(sub.isBestMac == 1, 1);
            if isempty(bi)
                bestMac = ''; nNoise = NaN; nFl = NaN; nFm = NaN; nFa = NaN;
            else
                bestMac = char(string(sub.channel(bi)));
                nNoise = getcol(sub,'nNoiseTrials',bi);
                nFl = getcol(sub,'nFinal_low',bi);
                nFm = getcol(sub,'nFinal_med',bi);
                nFa = getcol(sub,'nFinal_air',bi);
            end
            tot = sum([nFl nFm nFa], 'omitnan');
            if ~isfinite(tot) || tot == 0, status = 'excluded(allNoise)'; else, status = 'ok'; end
        end
        rows(i,:) = {char(T.subID(i)), id, T.sessNum(i), char(T.type(i)), char(T.group(i)), ...
                     status, nMac, bestMac, nNoise, nFl, nFm, nFa};
    end

    man = cell2table(rows, 'VariableNames', {'subID','sessID','sessNum','type','group', ...
        'status','nMacBP','bestMac','nNoiseTrials','nFinal_low','nFinal_med','nFinal_air'});
    writetable(man, manPath);
    fprintf('Wrote manifest %s (%d sessions; %d ok, %d noMacBP, %d excluded)\n', manPath, ...
        height(man), sum(strcmp(man.status,'ok')), sum(strcmp(man.status,'noMacBP')), ...
        sum(strcmp(man.status,'excluded(allNoise)')));
end

% ------------------------------------------------------------------
function v = local_getcol(sub, name, bi)
% return sub.(name)(bi) if the column exists, else NaN
    if ismember(name, sub.Properties.VariableNames)
        v = sub.(name)(bi);
    else
        v = NaN;
    end
end
