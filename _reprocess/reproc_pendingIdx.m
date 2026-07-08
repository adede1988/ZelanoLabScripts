function idx = reproc_pendingIdx(task)
% Curated session indices for TASK that are NOT yet in the done-marker. Used to
% set RUNSUBSET for the make step so it stays in sync with the main step's
% resume logic (never regenerate an intermediate for a session already finalized,
% which would clobber its good final).
    [cur, cfg] = reproc_curatedIdx(task);
    switch lower(char(task))
        case {'cuetask','cue'},              key = 'cue';
        case {'threshtask','thresh'},        key = 'thresh';
        case 'o15',                          key = 'O15';
        case {'breathingtask','breathing'},  key = 'breathing';
        otherwise, error('reproc_pendingIdx: unknown task %s', char(task));
    end
    doneIds = reproc_readDone(fullfile(reproc_root(), [key '_done.txt']));
    keep = true(size(cur));
    for i = 1:numel(cur)
        if any(strcmpi(cfg.sessionIDs{cur(i)}, doneIds)), keep(i) = false; end
    end
    idx = cur(keep);
end
