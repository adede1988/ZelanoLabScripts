function reproc_backup(task)
% Copy every curated final for TASK to E:\reprocBackup\<key>\ BEFORE the make
% step overwrites it. Skips a session already backed up (idempotent), so it is
% safe to re-run. NB: for O15 the runner backs up internally (no make step), so
% this is only needed for cue / thresh / breathing.
    [idx, cfg] = reproc_curatedIdx(task);
    [key, suffix] = taskFinal(task);
    dst = fullfile(reproc_root(), key);
    if ~isfolder(dst), mkdir(dst); end
    logf = fopen(fullfile(reproc_root(), [key '_backup.log']), 'a');
    logmsg(logf, sprintf('==== reproc_backup %s %s ====', key, datestr(now)));
    for s = idx(:)'
        id  = cfg.sessionIDs{s};
        src = fullfile(cfg.root{s}, id, 'preProc', [id suffix]);
        out = fullfile(dst, [id suffix]);
        try
            if exist(out, 'file'), logmsg(logf, sprintf('BKUP_EXISTS %s', id)); continue; end
            if exist(src, 'file')
                copyfile(src, out); logmsg(logf, sprintf('BKUP_OK %s', id));
            else
                logmsg(logf, sprintf('BKUP_NOFINAL %s', id));
            end
        catch ME
            logmsg(logf, sprintf('BKUP_FAIL %s : %s', id, ME.message));
        end
    end
    logmsg(logf, 'BACKUP_DONE'); if logf > 0, fclose(logf); end
end

function [key, suffix] = taskFinal(task)
    switch lower(char(task))
        case {'cuetask','cue'},              key = 'cue';       suffix = '_cueTaskPreproc.mat';
        case {'threshtask','thresh'},        key = 'thresh';    suffix = '_PEA_threshold_preproc.mat';
        case 'o15',                          key = 'O15';       suffix = '_O15preproc.mat';
        case {'breathingtask','breathing'},  key = 'breathing'; suffix = '_breathingPreproc.mat';
        otherwise, error('reproc_backup: unknown task %s', char(task));
    end
end

function logmsg(fid, msg)
    fprintf('%s\n', msg);
    if fid > 0, fprintf(fid, '%s\n', msg); end
end
