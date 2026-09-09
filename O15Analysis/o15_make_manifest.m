function o15_make_manifest()
% O15_MAKE_MANIFEST  Regenerate O15Task_data_manifest.csv from the current
%   FOOOF+QC summary CSV + the (fresh) session table, so the data-availability
%   table always reflects the live analysis state:
%     status = 'ok'                | macBP present, spectrograms made
%              'noMacBP'           | session has no macBP channels (skipped)
%              'excluded(noBase/allNoise)' | macBP present but no usable baseline
%                                            or no clean sniffs of any type
%   Derived (not hand-written) so it cannot drift from the rest of the report.

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo, 'O15Analysis'));
    groupDir = getenv('O15_GROUPDIR');
    if isempty(groupDir), groupDir = fullfile(labPaths().figPath, 'groupStatFigs'); end
    csvPath  = fullfile(groupDir, 'O15Task_fooof_summary.csv');
    manPath  = fullfile(groupDir, 'O15Task_data_manifest.csv');

    T = o15_session_table(false); T = T(T.onDisk & T.fresh, :);
    F = readtable(csvPath); Fsess = string(F.sessID);

    n = height(T);
    rows = cell(n, 13);
    for i = 1:n
        id = char(T.sessID(i));
        sel = strcmpi(Fsess, id);
        nMac = 0; bestMac = ''; nFs = NaN; nFf = NaN; nFc = NaN; baseFr = NaN;
        if ~any(sel)
            status = 'noMacBP';
        else
            sub = F(sel, :);
            nMac = height(sub);
            bi = find(sub.isBestMac == 1, 1);
            if isempty(bi)
                status = 'noMacBP';
            else
                bestMac = char(string(sub.channel(bi)));
                nFs = getcol(sub,'nFinal_start',bi); nFf = getcol(sub,'nFinal_free',bi);
                nFc = getcol(sub,'nFinal_confirm',bi); baseFr = getcol(sub,'baseFrames',bi);
                anyFinal = sum([nFs nFf nFc] > 0, 'omitnan') > 0;
                if ~(isfinite(baseFr) && baseFr > 0) || ~anyFinal
                    status = 'excluded(noBase/allNoise)';
                else
                    status = 'ok';
                end
            end
        end
        rows(i,:) = {char(T.subID(i)), id, T.sessNum(i), char(T.type(i)), char(T.group(i)), ...
                     status, nMac, bestMac, nFs, nFf, nFc, baseFr, char(T.modWhen(i))};
    end

    man = cell2table(rows, 'VariableNames', {'subID','sessID','sessNum','type','group', ...
        'status','nMacBP','bestMac','nFinal_start','nFinal_free','nFinal_confirm','baseFrames','modWhen'});
    writetable(man, manPath);
    fprintf('Wrote O15 manifest %s (%d sessions; %d ok, %d noMacBP, %d excluded)\n', manPath, ...
        height(man), sum(strcmp(man.status,'ok')), sum(strcmp(man.status,'noMacBP')), ...
        sum(contains(man.status,'excluded')));
end

function v = getcol(sub, name, bi)
    if ismember(name, sub.Properties.VariableNames), v = sub.(name)(bi); else, v = NaN; end
end
