function cue_make_manifest()
% CUE_MAKE_MANIFEST  Regenerate cueTask_data_manifest.csv from the current
%   FOOOF+QC summary CSV + the session table, so the data-availability table
%   always reflects the live analysis3 state:
%     status = 'ok'                | macBP present, spectrograms made
%              'noMacBP'           | session has no macBP channels (skipped)
%              'excluded(allNoise)'| macBP present but all trials flagged noisy
%   Derived (not hand-written) so it cannot drift from the rest of the report.

    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo); addpath(fullfile(repo, 'cueAnalysis'));
    groupDir = 'R:\Neurology\Zelano_Lab\Lab_Common\Adam\Dupi_processing\groupStatFigs';
    csvPath  = fullfile(groupDir, 'cueTask_fooof_summary.csv');
    manPath  = fullfile(groupDir, 'cueTask_data_manifest.csv');

    T = cue_session_table(false); T = T(T.onDisk, :);
    F = readtable(csvPath); Fsess = string(F.sessID);

    n = height(T);
    rows = cell(n, 11);
    for i = 1:n
        id = char(T.sessID(i));
        sel = strcmpi(Fsess, id);
        if ~any(sel)
            status = 'noMacBP'; nMac = 0; bestMac = ''; nNoise = NaN; nFt = NaN; nFf = NaN;
        else
            sub = F(sel, :);
            nMac = height(sub);
            bi = find(sub.isBestMac == 1, 1);
            if isempty(bi)
                bestMac = ''; nNoise = NaN; nFt = NaN; nFf = NaN;
            else
                bestMac = char(string(sub.channel(bi)));
                nNoise = sub.nNoiseTrials(bi); nFt = sub.nFinal_trialStart(bi); nFf = sub.nFinal_finalOnset(bi);
            end
            if ~isfinite(nFt) || nFt == 0, status = 'excluded(allNoise)'; else, status = 'ok'; end
        end
        rows(i,:) = {char(T.subID(i)), id, T.sessNum(i), char(T.type(i)), char(T.group(i)), ...
                     status, nMac, bestMac, nNoise, nFt, nFf};
    end

    man = cell2table(rows, 'VariableNames', {'subID','sessID','sessNum','type','group', ...
        'status','nMacBP','bestMac','nNoiseTrials','nFinal_trialStart','nFinal_finalOnset'});
    writetable(man, manPath);
    fprintf('Wrote manifest %s (%d sessions; %d ok, %d noMacBP, %d excluded)\n', manPath, ...
        height(man), sum(strcmp(man.status,'ok')), sum(strcmp(man.status,'noMacBP')), ...
        sum(strcmp(man.status,'excluded(allNoise)')));
end
