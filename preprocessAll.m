% preprocessAll.m
% ---------------------------------------------------------------------------
% Master driver: read dataTracking.xlsx, find every session whose raw data is
% extracted but that is NOT yet preprocessed (Data Preprocessed != 'X'), and run
% the appropriate task pipeline(s) for it.
%
% USAGE
%   preprocessAll                       % REPORT ONLY: print what's pending
%   setenv('PREPROCESS_RUN','1'); preprocessAll   % actually run the pending tasks
%
% A task's pipeline is its <task>_makeOutDat (if it has one) followed by its
% <task>PreProc_main. Those scripts already loop over all of their sessions and
% skip ones already done, so running them clears the remaining backlog for that
% task. They each begin with `clear`, so this driver recomputes its pending list
% fresh (via pendingReport) before each run rather than holding state across runs.
%
% NOTES
%  - GUESS sessions need interactive verification (paramCheck / onset gate) and
%    will halt a batch. Curate them first; pending guess rows are flagged below.
%  - Memory (16 GB): each task runs all its sessions in one MATLAB process. If you
%    hit "Out of memory", run one task at a time, or use the per-session
%    _dev/run_* harnesses (they clear big vars each iteration).
% ---------------------------------------------------------------------------

addpath(fileparts(mfilename('fullpath')));   % ensure labPaths/applyParams reachable

showReport(pendingReport());

if ~strcmp(getenv('PREPROCESS_RUN'), '1')
    fprintf(['\nReport only. To execute the pending pipelines:\n' ...
             '    setenv(''PREPROCESS_RUN'',''1''); preprocessAll\n']);
    return
end

% --- run pending tasks. Each `if` recomputes pendingReport() fresh because the
%     run() calls below execute `clear`. Order: makeOutDat (raw->intermediate)
%     then main (intermediate->final). ---
if pendingReport().breathingTask.n > 0
    run('breathingTask_makeOutDat.m');
    run('breathingTaskPreProc_main.m');
end
if pendingReport().cueTask.n > 0
    run('cueTask_makeOutDat.m');
    run('cueTaskPreProc_main.m');
end
if pendingReport().threshTask.n > 0
    run(fullfile('preproc', 'threshPreProc_makeOutDat.m'));
    run('threshPreProc_main.m');
end
if pendingReport().O15.n > 0     % O15 has no makeOutDat
    run('O15PreProc_main.m');
end

fprintf('\npreprocessAll: done. Re-run (report only) to confirm the backlog cleared.\n');


% =========================== local functions ===========================

function rep = pendingReport()
% Build, per task, the list of raw-extracted sessions not yet marked preprocessed.
    tasks = {'breathingTask', 'cueTask', 'threshTask', 'O15'};
    xlsx  = resolveSheet();
    C     = readcell(xlsx, 'Sheet', 'Sheet1');

    hRow = findHeaderRow(C);
    cSub = findCol(C(hRow,:), 'Subject ID');
    cTsk = findCol(C(hRow,:), 'Task');
    cPre = findCol(C(hRow,:), 'Data Preprocessed');

    rep = struct();
    for t = 1:numel(tasks)
        task = tasks{t};
        cfg  = applyParams(task, 'main');           % raw-extracted session list
        ids = {}; src = {};
        for s = 1:numel(cfg.sessionIDs)
            id  = cfg.sessionIDs{s};
            xval = lookupCell(C, hRow, cSub, cTsk, cPre, id, task);
            if ~strcmpi(strtrim(asChar(xval)), 'X')   % not yet preprocessed
                ids{end+1} = id; %#ok<AGROW>
                if numel(cfg.paramSource) >= s, src{end+1} = cfg.paramSource{s}; %#ok<AGROW>
                else, src{end+1} = ''; end %#ok<AGROW>
            end
        end
        rep.(task) = struct('ids', {ids}, 'src', {src}, 'n', numel(ids));
    end
end

function showReport(rep)
    tasks = fieldnames(rep);
    fprintf('\n==================== preprocessAll: pending ====================\n');
    total = 0;
    for t = 1:numel(tasks)
        r = rep.(tasks{t});
        nGuess = sum(strcmpi(r.src, 'guess'));
        fprintf('%-14s %2d pending  (%d guess -> need interactive verification)\n', ...
                tasks{t}, r.n, nGuess);
        for i = 1:r.n
            tag = ''; if strcmpi(r.src{i}, 'guess'), tag = '  [GUESS]'; end
            fprintf('                 - %s%s\n', r.ids{i}, tag);
        end
        total = total + r.n;
    end
    fprintf('----------------------------------------------------------------\n');
    fprintf('total pending: %d\n', total);
    fprintf('================================================================\n');
end

function xlsx = resolveSheet()
% Mirror applyParams: Admin master if it carries the param columns, else the
% repo-local copy next to labPaths.
    L = labPaths();
    xlsx = L.adminXlsx;
    if exist(xlsx, 'file') ~= 2 || ~headerHas(xlsx, 'datPre')
        local = fullfile(L.repo, 'dataTracking.xlsx');
        if exist(local, 'file') == 2, xlsx = local; end
    end
end

function tf = headerHas(xlsx, name)
    tf = false;
    try
        C = readcell(xlsx, 'Sheet', 'Sheet1');
        tf = findCol(C(findHeaderRow(C), :), name) > 0;
    catch
    end
end

function r = findHeaderRow(C)
    r = 0;
    for i = 1:size(C, 1)
        if findCol(C(i, :), 'Subject ID') > 0, r = i; return; end
    end
    if r == 0, error('preprocessAll:noHeader', 'No "Subject ID" header found.'); end
end

function idx = findCol(row, name)
    idx = 0;
    for c = 1:numel(row)
        v = row{c};
        if (ischar(v) || isstring(v)) && strcmpi(strtrim(char(v)), name), idx = c; return; end
    end
end

function v = lookupCell(C, hRow, cSub, cTsk, cPre, id, task)
% Data Preprocessed cell for the row matching this sessID + task (else blank).
    v = '';
    for r = hRow+1 : size(C, 1)
        sub = strtrim(asChar(C{r, cSub}));
        if strcmpi(sub, id) && strcmp(canonTask(C{r, cTsk}), taskKey(task))
            if cPre > 0 && cPre <= size(C, 2), v = C{r, cPre}; end
            return;
        end
    end
end

function s = asChar(v)
    if isa(v, 'missing'), s = ''; elseif ischar(v), s = v; elseif isstring(v), s = char(v);
    elseif isnumeric(v) || islogical(v), s = num2str(v); else, s = char(string(v)); end
end

function k = taskKey(task)
    s = lower(char(task)); s = s(~isspace(s));
    switch s
        case 'breathingtask', k = 'breathing';
        case 'cuetask',       k = 'cue';
        case 'threshtask',    k = 'thresh';
        case 'o15',           k = 'O15';
        otherwise,            k = '';
    end
end

function k = canonTask(t)
    if isa(t,'missing') || (isnumeric(t) && all(isnan(t(:)))) || isempty(t), k = ''; return; end
    s = lower(asChar(t)); s = s(~isspace(s));
    switch s
        case {'breathingtasks', 'wavebreathing', 'breathingtask'}, k = 'breathing';
        case 'odorcuetask', k = 'cue';
        case 'o15',         k = 'O15';
        case 'threshold',   k = 'thresh';
        otherwise,          k = '';
    end
end
