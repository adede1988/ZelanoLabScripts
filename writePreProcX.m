function writePreProcX(P, sessID, xlsxPath)
%WRITEPREPROCX  Mark a session's "Data Preprocessed" cell with an X.
%
%   writePreProcX(P, sessID)            % default dataTracking.xlsx
%   writePreProcX(P, sessID, xlsxPath)  % explicit spreadsheet path
%
%   Finds the row matching sessID (Subject ID) and P.task, writes 'X' into
%   that row's "Data Preprocessed" column, and saves the workbook in place.
%   Subject IDs repeat across tasks in the sheet, so the task is needed to
%   pick the right row.

    sheet = 'Sheet1';

    if nargin < 3 || isempty(xlsxPath)
        xlsxPath = resolveDefaultXlsx();
    end
    if exist(xlsxPath, 'file') ~= 2
        error('writePreProcX:noFile', 'Spreadsheet not found: %s', xlsxPath);
    end

    sessID = strtrim(char(string(sessID)));
    if isempty(sessID)
        error('writePreProcX:badSession', 'sessID must be a non-empty string.');
    end
    if ~isfield(P, 'task') || isempty(P.task)
        error('writePreProcX:noTask', 'P.task is required to identify the row.');
    end
    tkey = taskKey(P.task);

    % --- read the sheet, locate the header row (the one with 'Subject ID') ---
    C    = readcell(xlsxPath, 'Sheet', sheet);
    hRow = findRow(C, 'Subject ID');
    if hRow == 0
        error('writePreProcX:noHeader', 'No "Subject ID" header found in %s.', xlsxPath);
    end
    hdr = C(hRow, :);

    cSub  = findCol(hdr, 'Subject ID');
    cTask = findCol(hdr, 'Task');
    cPre  = findCol(hdr, 'Data Preprocessed');
    if cSub == 0 || cTask == 0 || cPre == 0
        error('writePreProcX:missingCols', ...
            'Sheet needs "Subject ID", "Task", and "Data Preprocessed" columns.');
    end

    % --- find the matching data row (Subject ID + task) ---
    dRow = 0;
    for r = hRow+1 : size(C, 1)
        sub = strtrim(char(string(C{r, cSub})));
        if strcmpi(sub, sessID) && strcmp(canonTask(C{r, cTask}), tkey)
            dRow = r; break;
        end
    end
    if dRow == 0
        error('writePreProcX:noSession', ...
            'No row found for session "%s" / task "%s".', sessID, char(string(P.task)));
    end

    % --- write the X and save ---
    ref = sprintf('%s%d', colLetter(cPre), dRow);
    writecell({'X'}, xlsxPath, 'Sheet', sheet, 'Range', ref);
    fprintf('writePreProcX: %s (%s) -> Data Preprocessed = X at %s\n', ...
            sessID, char(string(P.task)), ref);
end

% ============================ helpers ============================

function xlsxPath = resolveDefaultXlsx()
    L = labPaths();
    if exist(L.adminXlsx, 'file') == 2
        xlsxPath = L.adminXlsx;
    else
        xlsxPath = fullfile(L.repo, 'dataTracking.xlsx');
    end
end

function r = findRow(C, name)
    r = 0;
    for i = 1:size(C, 1)
        if findCol(C(i, :), name) > 0, r = i; return; end
    end
end

function idx = findCol(hdr, name)
    idx = 0;
    for c = 1:numel(hdr)
        v = hdr{c};
        if (ischar(v) || isstring(v)) && strcmpi(strtrim(char(v)), name)
            idx = c; return;
        end
    end
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

% --- task vocabulary: P.task and the sheet's Task column use different
%     names, so both map onto a common key before comparison. ---

function k = taskKey(task)   % P.task -> key
    s = lower(char(string(task)));
    s = s(~isspace(s));
    switch s
        case 'breathingtask', k = 'breathing';
        case 'cuetask',       k = 'cue';
        case 'threshtask',    k = 'thresh';
        case 'o15',           k = 'O15';
        otherwise,            k = '';
    end
end

function k = canonTask(t)    % sheet Task cell -> key
    if isBlank(t)
        k = ''; return;
    end
    s = lower(char(string(t)));
    s = s(~isspace(s));
    switch s
        case {'breathingtasks', 'wavebreathing', 'breathingtask'}
            k = 'breathing';
        case 'odorcuetask'
            k = 'cue';
        case 'o15'
            k = 'O15';
        case 'threshold'
            k = 'thresh';
        otherwise
            k = '';
    end
end

function tf = isBlank(v)
    if isa(v, 'missing'), tf = true; return; end
    if isnumeric(v) && (isempty(v) || all(isnan(v(:)))), tf = true; return; end
    if (ischar(v) || isstring(v)) && strlength(strtrim(string(v))) == 0, tf = true; return; end
    if iscell(v) && isempty(v), tf = true; return; end
    tf = false;
end
