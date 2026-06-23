function out = applyParams(task, sel, xlsxPath)
% applyParams  Single source of truth for session lists and per-session params.
%
%   Reads dataTracking.xlsx and returns either:
%     Mode A (cfg)  : a session-list config for a loop      (sel = 'makeOutDat'|'main')
%     Mode B (P)    : a per-session parameter struct        (sel = a session ID)
%
%   Replaces the hard-coded arrays in the *_makeOutDat / *_main scripts and the
%   entire parameter switch in getSessionParams_<task>.m.
%
%   Usage
%     cfg = applyParams(task, 'makeOutDat' | 'main' [, xlsxPath])   % Mode A
%     P   = applyParams(task, sessID                  [, xlsxPath]) % Mode B
%
%   task in {'breathingTask','cueTask','threshTask','O15'}
%   default xlsxPath = R:\Neurology\Zelano_Lab\Lab_Common\Admin\dataTracking.xlsx
%
%   Mode A cfg fields:
%     .sessionIDs (n x 1 cell)  .root (n x 1 cell)  .datPre (1 x k cell, fixed order)
%     .datPrei (1 x n)  .isNewStd (1 x n logical)  .newIDs (cell)
%     .rspIDX (1 x n)  .rspFlip (1 x n)  .paramSource (1 x n cell)
%
%   See CLAUDE.md / taskList.md for the full contract.

    if nargin < 3 || isempty(xlsxPath)
        xlsxPath = resolveDefaultXlsx();
    end

    % --- canonical roots (order is load-bearing: 1=Dupi 2=OBEControl 3=EEGbreathing) ---
    % Sourced from labPaths so they follow the machine's lab-common drive.
    L = labPaths();
    ROOT_DUPI = L.rootDupi;
    ROOT_OBE  = L.rootOBE;
    ROOT_EEG  = L.rootEEG;

    tkey = taskKey(task);
    if isempty(tkey)
        error('applyParams:badTask', 'Unknown task "%s".', char(string(task)));
    end

    C = readSheetCached(xlsxPath);

    % --- locate columns by trimmed header name (header row = 2) ---
    hdr = C(2, :);
    col = @(name) findCol(hdr, name);
    cSub  = col('Subject ID');
    cType = col('Type');
    cTask = col('Task');
    cRaw  = col('Raw Data Extracted');
    cPre  = col('datPre');
    cRsp  = col('rspIDX');
    cFlip = col('rspFlip');
    cHasE = col('hasEEG');
    cSpkC = col('spikeClean');
    cSpkT = col('spikeThresh');
    cSpkW = col('spikeWin');
    cMacR = col('macroRemove');
    cHasM = col('hasMacros');
    cRespT= col('respThresh');
    cBack = col('cuedBackBuff');
    cAdj  = col('adjWin');
    cBeat = col('beatSpec');
    cTtlR = col('ttlRemoveIdx');
    cTtlN = col('ttlNote');
    cNew  = col('isNewStd');
    cPS   = col('paramSource');

    % --- select rows: task matches AND raw extracted (not blank, not INCOMPLETE) ---
    data = C(3:end, :);
    nAll = size(data, 1);
    keep = false(nAll, 1);
    for r = 1:nAll
        if isBlank(data{r, cSub}), continue; end          % drop fully empty / no-ID rows
        if ~strcmp(canonTask(data{r, cTask}), tkey), continue; end
        rawv = data{r, cRaw};
        if isBlank(rawv), continue; end
        if strcmpi(strtrim(asChar(rawv)), 'INCOMPLETE'), continue; end
        keep(r) = true;
    end
    rows = data(keep, :);

    % --- dedupe by sessID (case-insensitive, keep first), stable sheet order ---
    n = size(rows, 1);
    seen = {};
    keep2 = false(n, 1);
    for r = 1:n
        id = strtrim(asChar(rows{r, cSub}));
        if any(strcmpi(id, seen)), continue; end
        seen{end+1} = id; %#ok<AGROW>
        keep2(r) = true;
    end
    rows = rows(keep2, :);
    n = size(rows, 1);

    selStr = strtrim(asChar(sel));
    isModeA = any(strcmpi(selStr, {'makeOutDat', 'main'}));

    if isModeA
        % ------------------------- Mode A: cfg -------------------------
        cfg = struct();
        cfg.sessionIDs = cell(n, 1);
        cfg.root       = cell(n, 1);
        cfg.rspIDX     = zeros(1, n);
        cfg.rspFlip    = zeros(1, n);
        cfg.isNewStd   = false(1, n);
        cfg.paramSource= cell(1, n);

        datPre = {ROOT_DUPI, ROOT_OBE, ROOT_EEG};
        for r = 1:n
            id = strtrim(asChar(rows{r, cSub}));
            cfg.sessionIDs{r} = id;

            rt = rootForRow(rows{r, cPre}, rows{r, cType}, L);
            cfg.root{r} = rt;

            cfg.rspIDX(r)  = num_or(rows{r, cRsp}, 1);
            cfg.rspFlip(r) = num_or(rows{r, cFlip}, 1);
            cfg.isNewStd(r)= bool_or(rows{r, cNew}, false);
            ps = rows{r, cPS};
            if isBlank(ps), cfg.paramSource{r} = ''; else, cfg.paramSource{r} = asChar(ps); end
        end

        % build datPre with any extra distinct roots appended after the fixed 3
        for r = 1:n
            if ~anyEqI(cfg.root{r}, datPre)
                datPre{end+1} = cfg.root{r}; %#ok<AGROW>
            end
        end
        cfg.datPre = datPre;

        cfg.datPrei = zeros(1, n);
        for r = 1:n
            cfg.datPrei(r) = idxEqI(cfg.root{r}, datPre);
        end

        cfg.newIDs = cfg.sessionIDs(cfg.isNewStd);

        out = cfg;
        return;
    end

    % ------------------------- Mode B: P -------------------------
    ri = 0;
    for r = 1:n
        if strcmpi(selStr, strtrim(asChar(rows{r, cSub})))
            ri = r; break;
        end
    end
    if ri == 0
        error('applyParams:noSession', ...
            'No %s row with raw-extracted data for session "%s".', task, selStr);
    end

    typeStudy = studyOf(rows{ri, cType});

    P = struct();
    P.task        = taskCallerKey(task);
    P.type        = typeStr(typeStudy);
    P.fs_target   = 500;
    P.debug       = false;
    P.computeResp = true;
    P.rspIDX      = num_or(rows{ri, cRsp}, 1);
    P.rspFlip     = num_or(rows{ri, cFlip}, 1);

    isO15 = strcmp(tkey, 'O15');
    P.hasEEG     = bool_or(rows{ri, cHasE}, ~isO15);   % default true except O15
    P.spikeClean = bool_or(rows{ri, cSpkC}, ~isO15);   % default true except O15
    P.spikeThresh= num_or(rows{ri, cSpkT}, 20);
    P.spikeWin   = num_or(rows{ri, cSpkW}, 11);
    P.macroRemove= list_or(rows{ri, cMacR}, []);
    P.paramSource= asChar(rows{ri, cPS});

    switch tkey
        case 'breathing'
            P.hasMacros = bool_or(rows{ri, cHasM}, true);
            beatSpec    = asChar(rows{ri, cBeat});
            if isBlank(rows{ri, cBeat}), beatSpec = '1,0,gt,3'; end
            P.beatSpec  = beatSpec;
            P.getBeats  = @(ECGz, beatSep) detectBeats(ECGz, beatSep, beatSpec);

        case 'cue'
            P.respThresh   = num_or(rows{ri, cRespT}, 500);
            P.cuedBackBuff = num_or(rows{ri, cBack}, 150);
            P.adjWin       = num_or(rows{ri, cAdj}, 500);
            P.ttlMap = struct( ...
                'cue',    {'cue','Cue','cueOnset'}, ...
                'target', {'targ','target','TargetOnset'}, ...
                'resp',   {'resp','response','button'} );

        case 'thresh'
            P.respThresh   = num_or(rows{ri, cRespT}, 500);
            P.cuedBackBuff = num_or(rows{ri, cBack}, 150);
            P.adjWin       = num_or(rows{ri, cAdj}, 500);
            P.ttlMap = struct('sniff', {'sniff'});

        case 'O15'
            P.respThresh   = num_or(rows{ri, cRespT}, 500);
            P.cuedBackBuff = num_or(rows{ri, cBack}, 150);
            P.adjWin       = num_or(rows{ri, cAdj}, 500);
            P.pd = struct('zthresh', -2, 'minPulseSamp', 200, ...
                          'maxPulseSamp', 1200, 'trialSplitSamp', 850);
            P.ttl = struct( ...
                'expectedTrialCount', 30, ...
                'removeTrialMarksIdx', list_or(rows{ri, cTtlR}, []), ...
                'note', stringOrEmpty(rows{ri, cTtlN}));
    end

    out = P;
end

% ======================= helpers =======================

function p = resolveDefaultXlsx()
% Prefer the lab Admin master IF it already carries the parameter columns;
% otherwise fall back to the repo-local param-enriched dataTracking.xlsx that
% lives next to applyParams.m. Resolution is cached for the session.
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
        p = adminP;   % last resort -> readSheetCached errors informatively
    end
    RESOLVED = p;
end

function tf = hasParamCols(p)
    tf = false;
    if exist(p, 'file') ~= 2, return; end
    try
        C = readcell(p, 'Sheet', 'Sheet1', 'Range', '2:2');   % header row only
        hdr = C(1, :);
        tf = any(cellfun(@(x) (ischar(x) || isstring(x)) && ...
                 strcmpi(strtrim(char(x)), 'datPre'), hdr));
    catch
        tf = false;
    end
end

function C = readSheetCached(xlsxPath)
    persistent CACHE
    d = dir(xlsxPath);
    if isempty(d)
        error('applyParams:noFile', 'Spreadsheet not found: %s', xlsxPath);
    end
    key = lower(strtrim(xlsxPath));
    if ~isempty(CACHE) && isfield(CACHE, 'key') && strcmp(CACHE.key, key) ...
            && CACHE.mtime == d.datenum
        C = CACHE.C;
        return;
    end
    C = readcell(xlsxPath, 'Sheet', 'Sheet1');
    CACHE = struct('key', key, 'mtime', d.datenum, 'C', {C});
end

function idx = findCol(hdr, name)
    idx = 0;
    for c = 1:numel(hdr)
        v = hdr{c};
        if ischar(v) || isstring(v)
            if strcmpi(strtrim(char(v)), name)
                idx = c; return;
            end
        end
    end
    if idx == 0
        error('applyParams:noColumn', 'Column "%s" not found in header row.', name);
    end
end

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

function s = stringOrEmpty(v)
    if isBlank(v), s = ""; else, s = string(asChar(v)); end
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
    % caller task key -> internal canon
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

function s = taskCallerKey(task)
    % normalize the caller's task spelling back to the canonical P.task value
    switch taskKey(task)
        case 'breathing', s = 'breathingTask';
        case 'cue',       s = 'cueTask';
        case 'thresh',    s = 'threshTask';
        case 'O15',       s = 'O15';
        otherwise,        s = asChar(task);
    end
end

function st = studyOf(typeVal)
    s = lower(asChar(typeVal));
    if contains(s, 'dupi')
        st = 'dupi';
    elseif contains(s, 'eeg')
        st = 'eeg';
    else
        st = 'obe';   % obecontrol / obe... / anything else
    end
end

function t = typeStr(study)
    switch study
        case 'dupi', t = 'Dupi';
        case 'eeg',  t = 'EEG';
        otherwise,   t = 'OBE';
    end
end

function rt = rootForRow(preVal, typeVal, L)
    if ~isBlank(preVal)
        rt = rebaseLabCommon(strtrim(asChar(preVal)), L);
        return;
    end
    switch studyOf(typeVal)
        case 'dupi', rt = L.rootDupi;
        case 'eeg',  rt = L.rootEEG;
        otherwise,   rt = L.rootOBE;
    end
end

function rt = rebaseLabCommon(rt, L)
% The sheet stores absolute datPre paths under the canonical Lab_Common prefix.
% Swap that prefix for this machine's labCommon so the sheet stays portable.
% No-op when they match (every machine that maps the share to R:).
    canon = L.labCommonCanon;
    if numel(rt) >= numel(canon) && strncmpi(rt, canon, numel(canon))
        rt = [L.labCommon, rt(numel(canon)+1:end)];
    end
end

function tf = anyEqI(s, lst)
    tf = false;
    for i = 1:numel(lst)
        if strcmpi(strtrim(s), strtrim(lst{i})), tf = true; return; end
    end
end

function idx = idxEqI(s, lst)
    idx = 0;
    for i = 1:numel(lst)
        if strcmpi(strtrim(s), strtrim(lst{i})), idx = i; return; end
    end
end
