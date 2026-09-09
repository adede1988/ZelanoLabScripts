function info = writeParams(P, sessID, xlsxPath, varargin)
%WRITEPARAMS  Persist an edited per-session P struct back to dataTracking.xlsx.
%
%   The inverse of applyParams(task, sessID) (Mode B). While you are
%   processing a session you may tweak parameters in P; calling
%
%       writeParams(P, S.id)
%
%   writes those values into the matching row of dataTracking.xlsx so the
%   spreadsheet stays the single source of truth. Only the parameter columns
%   that applyParams reads are touched; identity columns (Subject ID, Type,
%   Task, Raw Data Extracted) are never modified, and only cells whose value
%   actually changed are written.
%
%   Usage
%     writeParams(P, sessID)
%     writeParams(P, sessID, xlsxPath)
%     info = writeParams(P, sessID, xlsxPath, 'Name', value, ...)
%
%   The task is taken from P.task (e.g. 'cueTask'); override with 'Task'.
%   Default xlsxPath mirrors applyParams: the lab Admin master if it already
%   carries the param columns, else the repo-local dataTracking.xlsx next to
%   this file.
%
%   Name-value options
%     'Verbose'          (true)     print a one-line summary of each change
%     'DryRun'           (false)    compute & print changes but do not write
%     'Sheet'            ('Sheet1') worksheet name
%     'Task'             ('')       task override if P.task is absent/wrong
%     'CreateMissingCols'(false)    append a new header column for any param
%                                   column not present in the sheet
%
%   Returns (optional) info struct:
%     .xlsxPath .sheet .row .session .task
%     .written            (cell of structs: col, ref, old, new)
%     .unchanged          (cell of column names left as-is)
%     .skippedMissingCol  (cell of column names not in the sheet)
%     .createdCol         (cell of column names newly appended)
%
%   Parameter columns written (when present in P and in the sheet):
%     rspIDX rspFlip hasEEG spikeClean spikeThresh spikeWin macroRemove
%     paramSource hasMacros beatSpec respThresh cuedBackBuff adjWin
%     ttlRemoveIdx (<- P.ttl.removeTrialMarksIdx) ttlNote (<- P.ttl.note)
%   Non-serialisable fields (getBeats, ttlMap, pd) are ignored by design.

    % ---------------- options ----------------
    opt.Verbose           = true;
    opt.DryRun            = false;
    opt.Sheet             = 'Sheet1';
    opt.Task              = '';
    opt.CreateMissingCols = false;
    opt = parseOpts(opt, varargin);
    sheet = opt.Sheet;

    if nargin < 3 || isempty(xlsxPath)
        xlsxPath = resolveDefaultXlsx();
    end

    sessID = strtrim(asChar(sessID));
    if isempty(sessID)
        error('writeParams:badSession', 'sessID must be a non-empty string.');
    end

    % ---------------- task ----------------
    if ~isempty(opt.Task)
        taskName = asChar(opt.Task);
    elseif isfield(P, 'task') && ~isBlank(P.task)
        taskName = asChar(P.task);
    else
        error('writeParams:noTask', ...
            'Cannot determine task: P has no .task field and no ''Task'' option given.');
    end
    tkey = taskKey(taskName);
    if isempty(tkey)
        error('writeParams:badTask', 'Unknown task "%s".', taskName);
    end

    % ---------------- read sheet ----------------
    if exist(xlsxPath, 'file') ~= 2
        error('writeParams:noFile', 'Spreadsheet not found: %s', xlsxPath);
    end
    C = readcell(xlsxPath, 'Sheet', sheet);

    % ---------------- locate header row (row containing 'Subject ID') ----------------
    hRow = 0;
    for r = 1:size(C, 1)
        if rowHasHeader(C(r, :), 'Subject ID')
            hRow = r; break;
        end
    end
    if hRow == 0
        error('writeParams:noHeader', ...
            'Could not find a header row containing "Subject ID" in %s.', xlsxPath);
    end
    hdr = C(hRow, :);

    cSub  = findColSoft(hdr, 'Subject ID');
    cTask = findColSoft(hdr, 'Task');
    cRaw  = findColSoft(hdr, 'Raw Data Extracted');
    if cSub == 0 || cTask == 0
        error('writeParams:missingKeyCols', ...
            'Sheet must contain "Subject ID" and "Task" columns.');
    end

    % ---------------- find target row (mirror applyParams selection) ----------------
    dRow = 0;
    for r = hRow+1 : size(C, 1)
        sub = strtrim(asChar(C{r, cSub}));
        if isempty(sub), continue; end
        if ~strcmp(canonTask(C{r, cTask}), tkey), continue; end
        if cRaw > 0
            rawv = C{r, cRaw};
            if isBlank(rawv), continue; end
            if strcmpi(strtrim(asChar(rawv)), 'INCOMPLETE'), continue; end
        end
        if strcmpi(sub, sessID)
            dRow = r; break;
        end
    end
    if dRow == 0
        error('writeParams:noSession', ...
            ['No %s row with raw-extracted data found for session "%s". ' ...
             'writeParams only updates rows that applyParams would read.'], ...
             taskName, sessID);
    end

    % ---------------- safety: confirm coordinate alignment ----------------
    % readcell maps C(i,j) -> Excel cell (i,j) only when the used range starts
    % at A1. Re-read the Subject ID cell at the computed coordinate and confirm
    % it matches before writing anything, so a shifted layout can't corrupt data.
    subRef = cellRef(cSub, dRow);
    try
        chk    = readcell(xlsxPath, 'Sheet', sheet, 'Range', subRef);
        chkVal = strtrim(asChar(chk{1}));
    catch ME
        error('writeParams:verifyFailed', ...
            'Could not re-read %s for verification: %s', subRef, ME.message);
    end
    if ~strcmpi(chkVal, sessID)
        error('writeParams:misaligned', ...
            ['Cell %s reads "%s" but expected "%s". The sheet does not map 1:1 ' ...
             'to cell coordinates (leading blank rows/columns?). Aborting to ' ...
             'avoid corrupting data.'], subRef, chkVal, sessID);
    end

    % ---------------- build writable entries from P ----------------
    % {colName, kind, valueGetter}; getter is only evaluated when present.
    E = {};
    E = addEntry(E, isfield(P,'rspIDX'),      'rspIDX',       'num',  @() P.rspIDX);
    E = addEntry(E, isfield(P,'rspFlip'),     'rspFlip',      'num',  @() P.rspFlip);
    E = addEntry(E, isfield(P,'hasEEG'),      'hasEEG',       'bool', @() P.hasEEG);
    E = addEntry(E, isfield(P,'spikeClean'),  'spikeClean',   'bool', @() P.spikeClean);
    E = addEntry(E, isfield(P,'spikeThresh'), 'spikeThresh',  'num',  @() P.spikeThresh);
    E = addEntry(E, isfield(P,'spikeWin'),    'spikeWin',     'num',  @() P.spikeWin);
    E = addEntry(E, isfield(P,'macroRemove'), 'macroRemove',  'list', @() P.macroRemove);
    E = addEntry(E, isfield(P,'paramSource'), 'paramSource',  'str',  @() P.paramSource);
    E = addEntry(E, isfield(P,'hasMacros'),   'hasMacros',    'bool', @() P.hasMacros);
    E = addEntry(E, isfield(P,'beatSpec'),    'beatSpec',     'str',  @() P.beatSpec);
    E = addEntry(E, isfield(P,'respThresh'),  'respThresh',   'num',  @() P.respThresh);
    E = addEntry(E, isfield(P,'cuedBackBuff'),'cuedBackBuff', 'num',  @() P.cuedBackBuff);
    E = addEntry(E, isfield(P,'adjWin'),      'adjWin',       'num',  @() P.adjWin);
    if isfield(P, 'ttl') && isstruct(P.ttl)
        E = addEntry(E, isfield(P.ttl,'removeTrialMarksIdx'), ...
                        'ttlRemoveIdx', 'list', @() P.ttl.removeTrialMarksIdx);
        E = addEntry(E, isfield(P.ttl,'note'), ...
                        'ttlNote',      'str',  @() P.ttl.note);
    end

    % ---------------- diff & write ----------------
    info = struct('xlsxPath', xlsxPath, 'sheet', sheet, 'row', dRow, ...
                  'session', sessID, 'task', taskName, ...
                  'written', {{}}, 'unchanged', {{}}, ...
                  'skippedMissingCol', {{}}, 'createdCol', {{}});

    lastCol    = lastNonBlankCol(hdr);   % running tail for CreateMissingCols
    createdCol = {};                     % {name, colIdx} resolved this call

    for k = 1:size(E, 1)
        colName = E{k, 1};
        kind    = E{k, 2};
        val     = E{k, 3};

        col = findColSoft(hdr, colName);
        if col == 0
            col = lookupCreated(createdCol, colName);
        end

        if col == 0
            if opt.CreateMissingCols
                lastCol = lastCol + 1;
                col     = lastCol;
                if ~opt.DryRun
                    writeOneCell(xlsxPath, sheet, cellRef(col, hRow), colName);
                end
                createdCol(end+1, :) = {colName, col};        %#ok<AGROW>
                info.createdCol{end+1} = colName;             %#ok<AGROW>
                if opt.Verbose
                    fprintf('%s+ new column "%s" at %s\n', dryTag(opt), ...
                            colName, colLetter(col));
                end
            else
                info.skippedMissingCol{end+1} = colName;      %#ok<AGROW>
                if opt.Verbose
                    warning('writeParams:noColumn', ...
                        'Column "%s" not in sheet; skipping (use CreateMissingCols to add).', ...
                        colName);
                end
                continue;
            end
        end

        cur = [];
        if col <= size(C, 2) && dRow <= size(C, 1)
            cur = C{dRow, col};
        end

        [newDisp, isChanged, oldDisp] = diffCell(cur, kind, val);
        if ~isChanged
            info.unchanged{end+1} = colName;                  %#ok<AGROW>
            continue;
        end

        ref = cellRef(col, dRow);
        if ~opt.DryRun
            writeOneCell(xlsxPath, sheet, ref, encodeCell(kind, val));
        end
        info.written{end+1} = struct('col', colName, 'ref', ref, ...
                                     'old', oldDisp, 'new', newDisp);   %#ok<AGROW>
        if opt.Verbose
            fprintf('%s%s/%s  %-13s %s -> %s\n', dryTag(opt), ...
                    taskName, sessID, colName, oldDisp, newDisp);
        end
    end

    if opt.Verbose && isempty(info.written) && isempty(info.createdCol)
        fprintf('writeParams: %s/%s  no changes.\n', taskName, sessID);
    end
end

% ============================ helpers ============================

function opt = parseOpts(opt, args)
    if mod(numel(args), 2) ~= 0
        error('writeParams:badOpts', 'Options must be name-value pairs.');
    end
    fn = fieldnames(opt);
    for i = 1:2:numel(args)
        name = args{i};
        if ~(ischar(name) || isstring(name))
            error('writeParams:badOptName', 'Option names must be strings.');
        end
        m = find(strcmpi(char(name), fn), 1);
        if isempty(m)
            error('writeParams:unknownOpt', 'Unknown option "%s".', char(name));
        end
        opt.(fn{m}) = args{i+1};
    end
end

function t = dryTag(opt)
    if opt.DryRun, t = '[dry-run] '; else, t = ''; end
end

function p = resolveDefaultXlsx()
% Mirror applyParams: prefer the lab Admin master only if it carries the
% parameter columns; otherwise use the repo-local param-enriched copy next
% to this file. Cached for the session.
    persistent RESOLVED
    if ~isempty(RESOLVED), p = RESOLVED; return; end
    L = labPaths();
    adminP = L.adminXlsx;
    localP = fullfile(L.repo, 'dataTracking.xlsx');
    if hasParamCols(adminP)
        p = adminP;
    elseif exist(localP, 'file') == 2
        p = localP;
    else
        p = adminP;   % last resort -> errors informatively on read
    end
    RESOLVED = p;
end

function tf = hasParamCols(p)
    tf = false;
    if exist(p, 'file') ~= 2, return; end
    try
        C   = readcell(p, 'Sheet', 'Sheet1');
        for r = 1:size(C, 1)
            if rowHasHeader(C(r, :), 'datPre'), tf = true; return; end
        end
    catch
        tf = false;
    end
end

function tf = rowHasHeader(rowCells, name)
    tf = false;
    for c = 1:numel(rowCells)
        v = rowCells{c};
        if (ischar(v) || isstring(v)) && strcmpi(strtrim(char(v)), name)
            tf = true; return;
        end
    end
end

function idx = findColSoft(hdr, name)
    idx = 0;
    for c = 1:numel(hdr)
        v = hdr{c};
        if (ischar(v) || isstring(v)) && strcmpi(strtrim(char(v)), name)
            idx = c; return;
        end
    end
end

function n = lastNonBlankCol(hdr)
    n = 0;
    for c = 1:numel(hdr)
        if ~isBlank(hdr{c}), n = c; end
    end
end

function col = lookupCreated(createdCol, name)
    col = 0;
    for i = 1:size(createdCol, 1)
        if strcmpi(createdCol{i, 1}, name), col = createdCol{i, 2}; return; end
    end
end

function E = addEntry(E, present, colName, kind, valFcn)
    if present
        E(end+1, :) = {colName, kind, valFcn()};
    end
end

function ref = cellRef(col, row)
    ref = sprintf('%s%d', colLetter(col), row);
end

function s = colLetter(n)
    s = '';
    while n > 0
        r = mod(n - 1, 26);
        s = [char('A' + r) s]; %#ok<AGROW>
        n = floor((n - 1) / 26);
    end
    if isempty(s), s = 'A'; end
end

function [newDisp, changed, oldDisp] = diffCell(cur, kind, val)
% Compare new value against the current cell using the same read semantics
% as applyParams, so formatting-only differences (e.g. 500 vs "500") are not
% counted as changes.
    switch kind
        case 'num'
            newV    = double(val);
            newDisp = num2str(newV);
            if isBlank(cur), oldDisp = '<blank>'; changed = true; return; end
            oldV    = num_or(cur, NaN);
            oldDisp = num2str(oldV);
            changed = ~(isfinite(oldV) && abs(oldV - newV) <= 1e-9);

        case 'bool'
            newV    = logical(val);
            newDisp = boolStr(newV);
            if isBlank(cur), oldDisp = '<blank>'; changed = true; return; end
            oldV    = bool_or(cur, ~newV);
            oldDisp = boolStr(oldV);
            changed = (oldV ~= newV);

        case 'list'
            newV    = double(val(:)).';
            newDisp = ['[' serializeList(newV) ']'];
            if isBlank(cur), oldDisp = '<blank>'; changed = ~isempty(newV); return; end
            oldV    = list_or(cur, []);
            oldDisp = ['[' serializeList(oldV) ']'];
            changed = ~isequal(oldV(:).', newV);

        case 'str'
            newV    = strtrim(asChar(val));
            newDisp = ['"' newV '"'];
            if isBlank(cur), oldS = ''; else, oldS = strtrim(asChar(cur)); end
            oldDisp = ['"' oldS '"'];
            changed = ~strcmp(oldS, newV);

        otherwise
            error('writeParams:badKind', 'Unknown kind "%s".', kind);
    end
end

function v = encodeCell(kind, val)
% Convert a P value into what gets written into the cell.
    switch kind
        case 'num',  v = double(val);
        case 'bool', v = logical(val);
        case 'list', v = serializeList(double(val(:)).');   % '' for empty -> blank-ish
        case 'str',  v = strtrim(asChar(val));
        otherwise,   error('writeParams:badKind', 'Unknown kind "%s".', kind);
    end
end

function writeOneCell(xlsxPath, sheet, ref, val)
% Update a single cell in place. writecell with a Range preserves the rest of
% the workbook (other cells and sheets). Falls back if AutoFitWidth is
% unsupported on older MATLAB.
    try
        writecell({val}, xlsxPath, 'Sheet', sheet, 'Range', ref, 'AutoFitWidth', false);
    catch e1 %#ok<NASGU>
        try
            writecell({val}, xlsxPath, 'Sheet', sheet, 'Range', ref);
        catch e2
            error('writeParams:writeFailed', ...
                ['Failed to write %s in %s.\nMake sure the file is closed in ' ...
                 'Excel and you have write access (network drive). Error: %s'], ...
                 ref, xlsxPath, e2.message);
        end
    end
end

function s = serializeList(v)
    v = v(:).';
    if isempty(v), s = ''; return; end
    parts = arrayfun(@(x) strtrim(num2str(x)), v, 'UniformOutput', false);
    s = strjoin(parts, ',');
end

function s = boolStr(b)
    if b, s = 'TRUE'; else, s = 'FALSE'; end
end

% ---- shared read primitives (kept identical to applyParams) ----

function tf = isBlank(v)
    if isa(v, 'missing'), tf = true; return; end
    if isnumeric(v) && isempty(v), tf = true; return; end
    if isnumeric(v) && isscalar(v) && isnan(v), tf = true; return; end
    if (ischar(v) || isstring(v)) && strlength(strtrim(string(v))) == 0, tf = true; return; end
    if iscell(v) && isempty(v), tf = true; return; end
    tf = false;
end

function s = asChar(v)
    if isa(v, 'missing'), s = ''; return; end
    if ischar(v), s = v; return; end
    if isstring(v), s = char(v); return; end
    if isnumeric(v) || islogical(v)
        if isscalar(v), s = num2str(v); else, s = mat2str(v); end
        return;
    end
    s = char(string(v));
end

function x = num_or(v, d)
    if isBlank(v), x = d; return; end
    if isnumeric(v) || islogical(v), x = double(v); return; end
    x = str2double(asChar(v));
    if isnan(x), x = d; end
end

function tf = bool_or(v, d)
    if isBlank(v), tf = d; return; end
    if islogical(v), tf = logical(v); return; end
    if isnumeric(v), tf = (v ~= 0); return; end
    s = lower(strtrim(asChar(v)));
    tf = ismember(s, {'true', '1', 'yes'});
end

function r = list_or(v, d)
    if isBlank(v), r = d; return; end
    if isnumeric(v), r = double(v(:))'; return; end
    if islogical(v), r = double(v(:))'; return; end
    s = asChar(v);
    try
        r = str2num(['[' s ']']); %#ok<ST2NM>
    catch
        r = d; return;
    end
    if isempty(r), r = d; end
end

function k = canonTask(t)
    if isBlank(t), k = ''; return; end
    s = lower(asChar(t));
    s = s(~isspace(s));
    switch s
        case {'breathingtasks', 'wavebreathing', 'breathingtask'}
            k = 'breathing';
        case {'odorcuetask'}
            k = 'cue';
        case {'o15'}
            k = 'O15';
        case {'threshold'}
            k = 'thresh';
        otherwise
            k = '';
    end
end

function k = taskKey(task)
    s = lower(asChar(task));
    s = s(~isspace(s));
    switch s
        case 'breathingtask', k = 'breathing';
        case 'cuetask',       k = 'cue';
        case 'threshtask',    k = 'thresh';
        case 'o15',           k = 'O15';
        otherwise,            k = '';
    end
end
