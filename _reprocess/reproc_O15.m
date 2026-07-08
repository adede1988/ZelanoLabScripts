function reproc_O15(idxList)
% REPROC_O15  Portable forced re-process of O15 sessions from raw, overwriting
% the R: finals. Mirrors O15PreProc_main.m / run_O15.m but: uses labPaths (runs
% on any machine), backs up each existing final before overwrite, processes only
% CURATED sessions (paramSource ~= 'guess'), isolates failures per session, and
% logs a SUCCESS/FAIL summary (incl. a manOnset sanity check).
%
%   reproc_O15()          % all curated O15 sessions
%   reproc_O15(idxList)   % only these session indices (curated ones among them)
%
% O15 has no makeOutDat step: assemble_outDat_all loads raw_O15 + detect_ttls_O15
% directly, so this single pass IS "from raw".

    L = labPaths();
    addpath(genpath(L.repo));
    addpath(genpath(L.eeglab));
    figPath = L.figPath;
    EEGLOC  = readtable(L.eegLocCsv);
    set(0, 'defaultfigurevisible', 'off');

    backupRoot = fullfile(reprocRoot(), 'O15');
    if ~isfolder(backupRoot), mkdir(backupRoot); end
    logf = fopen(fullfile(reprocRoot(), 'O15_run.log'), 'a');
    logmsg(logf, sprintf('==== reproc_O15 start %s ====', datestr(now)));

    donePath = fullfile(reprocRoot(), 'O15_done.txt');
    doneIds  = readDone(donePath);

    cfg = applyParams('O15', 'main');
    if nargin < 1 || isempty(idxList), idxList = 1:numel(cfg.sessionIDs); end

    for s = idxList(:)'
        id = cfg.sessionIDs{s};
        if strcmpi(cfg.paramSource{s}, 'guess')
            logmsg(logf, sprintf('GUESS_SKIP_O15 %s', id)); continue;
        end
        if any(strcmpi(id, doneIds))
            logmsg(logf, sprintf('ALREADY_DONE_O15 %s (skip; delete O15_done.txt to force)', id)); continue;
        end
        try
            S = struct('id', id, 'root', cfg.root{s}, 'figPath', figPath);
            finalMat = fullfile(S.root, S.id, 'preProc', [S.id '_O15preproc.mat']);
            if exist(finalMat, 'file')
                copyfile(finalMat, fullfile(backupRoot, [S.id '_O15preproc.mat']));
            end
            if ~isfolder(fullfile(S.figPath, S.id)), mkdir(fullfile(S.figPath, S.id)); end

            P = applyParams('O15', S.id);
            [outDat, raw, TTL] = assemble_outDat_all(S, P);
            outDat.TTL = TTL;
            outDat = downsample_data(outDat, P.fs_target);
            if P.hasEEG, outDat = preprocess_eeg(outDat, EEGLOC, P); end
            outDat = preprocess_macros(outDat, P);
            outDat.moreThan1 = 1;
            outDat.rspIDX = P.rspIDX; outDat.rspFlip = P.rspFlip;

            R = preprocess_respiration_wholetrace(outDat);
            sniffs = detect_sniffs_from_TTLs(R, P, outDat);
            outDat.behDat = build_behavior_table_O15(sniffs, raw.beh);
            outDat = refine_onsets_with_phase(outDat, R, P);
            plot_sniff_epochs(outDat, R);

            if ~isfolder(fullfile(outDat.OGdataDir, 'preProc')), mkdir(fullfile(outDat.OGdataDir, 'preProc')); end
            save(fullfile(outDat.OGdataDir, 'preProc', [outDat.sessID '_O15preproc.mat']), 'outDat', '-v7.3');

            hasMan = ismember('manOnset', outDat.behDat.Properties.VariableNames);
            logmsg(logf, sprintf('SUCCESS_O15 %s (manOnset=%d nBeh=%d nCh=%d fs=%d)', ...
                id, hasMan, height(outDat.behDat), numel(outDat.labels), outDat.fs));
            markDone(donePath, id);
        catch ME
            logmsg(logf, sprintf('FAIL_O15 %s : %s', id, ME.message));
        end
        clear outDat raw R sniffs TTL P
    end
    logmsg(logf, 'RUN_O15_DONE'); fclose(logf);
end

function r = reprocRoot()
% Backup/log root: E:\reprocBackup on the lab desktop, else <repo>\..\reprocBackup.
    if isfolder('E:\'), r = 'E:\reprocBackup';
    else, L = labPaths(); r = fullfile(fileparts(L.repo), 'reprocBackup'); end
    if ~isfolder(r), mkdir(r); end
end

function logmsg(fid, msg)
    fprintf('%s\n', msg);
    if fid > 0, fprintf(fid, '%s\n', msg); end
end

function ids = readDone(p)
    ids = {};
    if exist(p, 'file') == 2
        ids = strtrim(string(splitlines(fileread(p))));
        ids = cellstr(ids(strlength(ids) > 0));
    end
end

function markDone(p, id)
    fid = fopen(p, 'a'); if fid > 0, fprintf(fid, '%s\n', id); fclose(fid); end
end
