function dumpPreprocStructs()
% DUMPPREPROCSTRUCTS  Inspect the FINAL preprocessed .mat structures for each task.
%   Writes a human-readable description of the saved struct (fields, classes,
%   sizes, sample values, behDat table schema) to _dev/struct_dump.txt.
%
%   Memory-safe: enumerates fields for EVERY session cheaply via h5info (no load),
%   then deep-loads exactly ONE "rich" session per task (hasEEG, +hasMacros for
%   breathing) for value samples + table content, clearing it before the next task.

    here    = fileparts(mfilename('fullpath'));
    repo    = fileparts(here);
    addpath(repo); addpath(here);

    outTxt  = fullfile(here, 'struct_dump.txt');
    fid     = fopen(outTxt, 'w');
    cleaner = onCleanup(@() fclose(fid));
    P2 = @(varargin) fprintf(fid, varargin{:});   % to file
    say = @(varargin) fprintf(1, varargin{:});    % to console

    tasks = { ...
        'breathingTask', @(id) ['_breathingPreproc.mat']; ...
        'cueTask',       @(id) ['_cueTaskPreproc.mat']; ...
        'threshTask',    @(id) ['_PEA_threshold_preproc.mat']; ...
        'O15',           @(id) ['_O15preproc.mat'] };

    for ti = 1:size(tasks,1)
        task = tasks{ti,1};
        sfx  = tasks{ti,2}('');
        P2('\n\n########################################################\n');
        P2('# TASK: %s   (final file suffix: <id>%s)\n', task, sfx);
        P2('########################################################\n');
        say('=== %s ===\n', task);

        try
            cfg = applyParams(task, 'main');
        catch ME
            P2('  applyParams failed: %s\n', ME.message); continue;
        end
        ids   = cfg.sessionIDs;
        roots = cfg.root;
        P2('  %d sessions listed for this task.\n', numel(ids));

        % ---- cheap field-presence sweep over all sessions (h5info) ----
        allFields = {};
        fieldCount = containers.Map('KeyType','char','ValueType','double');
        existing = {};   % {id, root, fullpath, varname}
        for s = 1:numel(ids)
            id = ids{s}; root = roots{s};
            f  = findFinal(root, id, sfx);
            if isempty(f), continue; end
            try
                vn = topVar(f);
                flds = h5FieldNames(f, vn);
            catch
                flds = {};
            end
            existing(end+1,:) = {id, root, f, ''}; %#ok<AGROW>
            for k = 1:numel(flds)
                if ~ismember(flds{k}, allFields), allFields{end+1} = flds{k}; end %#ok<AGROW>
                if isKey(fieldCount, flds{k}), fieldCount(flds{k}) = fieldCount(flds{k})+1;
                else, fieldCount(flds{k}) = 1; end
            end
        end
        P2('  %d sessions have a final file on disk.\n', size(existing,1));
        if ~isempty(allFields)
            P2('\n  -- field presence across sessions (field : #sessions with it) --\n');
            ks = sort(allFields);
            for k = 1:numel(ks)
                P2('     %-22s : %d\n', ks{k}, fieldCount(ks{k}));
            end
        end
        if isempty(existing)
            P2('  (no final files found on disk for %s)\n', task); continue;
        end

        % ---- choose a RICH session to deep-dump ----
        pick = chooseRich(task, existing);
        P2('\n  -- DEEP DUMP of rich session: %s --\n', existing{pick,1});
        say('  deep-dumping %s ...\n', existing{pick,1});
        f  = existing{pick,3};
        vn = topVar(f);
        try
            Sv = load(f);
            od = Sv.(vn);
            clear Sv;
            describe(fid, od, vn, 0, task);
            clear od;
        catch ME
            P2('  DEEP LOAD FAILED (%s): %s\n', existing{pick,1}, ME.message);
            P2('  Falling back to h5info structure listing:\n');
            dumpH5(fid, f, vn);
        end
    end
    say('\nDONE. Wrote %s\n', outTxt);
end

% ======================================================================
function f = findFinal(root, id, sfx)
    f = '';
    d = fullfile(root, id, 'preProc');
    if ~isfolder(d), return; end
    cand = fullfile(d, [id sfx]);
    if exist(cand, 'file'), f = cand; return; end
    % case-insensitive / fuzzy fallback
    lst = dir(fullfile(d, [id '*.mat']));
    for i = 1:numel(lst)
        if strcmpi(lst(i).name, [id sfx]), f = fullfile(d, lst(i).name); return; end
    end
end

function vn = topVar(f)
    w = who('-file', f);
    pref = {'outDat','chanDat','out'};
    vn = '';
    for i = 1:numel(pref)
        if any(strcmp(w, pref{i})), vn = pref{i}; return; end
    end
    if ~isempty(w), vn = w{1}; end
end

function flds = h5FieldNames(f, vn)
    info = h5info(f);
    flds = {};
    g = findGroup(info, ['/' vn]);
    if isempty(g), return; end
    if isfield(g,'Groups') && ~isempty(g.Groups)
        for i = 1:numel(g.Groups)
            nm = g.Groups(i).Name; flds{end+1} = lastTok(nm); %#ok<AGROW>
        end
    end
    if isfield(g,'Datasets') && ~isempty(g.Datasets)
        for i = 1:numel(g.Datasets)
            flds{end+1} = g.Datasets(i).Name; %#ok<AGROW>
        end
    end
end

function g = findGroup(info, name)
    g = [];
    if strcmp(info.Name, name), g = info; return; end
    if isfield(info,'Groups')
        for i = 1:numel(info.Groups)
            g = findGroup(info.Groups(i), name);
            if ~isempty(g), return; end
        end
    end
end

function t = lastTok(nm)
    parts = strsplit(nm, '/');
    t = parts{end};
end

function pick = chooseRich(task, existing)
    pick = 1; best = -1;
    for i = 1:size(existing,1)
        id = existing{i,1};
        sc = 0;
        try
            Pp = applyParams(task, id);
            if isfield(Pp,'hasEEG') && Pp.hasEEG, sc = sc + 2; end
            if isfield(Pp,'hasMacros') && Pp.hasMacros, sc = sc + 1; end
        catch
        end
        if sc > best, best = sc; pick = i; end
        if best >= 3, break; end   % good enough
    end
end

function dumpH5(fid, f, vn)
    info = h5info(f);
    g = findGroup(info, ['/' vn]);
    if isempty(g), fprintf(fid,'   (no group /%s)\n', vn); return; end
    walkH5(fid, g, 1);
end

function walkH5(fid, g, depth)
    pad = repmat('  ',1,depth);
    if isfield(g,'Datasets')
        for i = 1:numel(g.Datasets)
            ds = g.Datasets(i);
            sz = '?'; try, sz = mat2str(ds.Dataspace.Size); catch, end
            fprintf(fid, '%s%-20s [%s] %s\n', pad, ds.Name, sz, ds.Datatype.Class);
        end
    end
    if isfield(g,'Groups')
        for i = 1:numel(g.Groups)
            fprintf(fid, '%s%s/\n', pad, lastTok(g.Groups(i).Name));
            if depth < 3, walkH5(fid, g.Groups(i), depth+1); end
        end
    end
end

% ======================================================================
function describe(fid, v, name, depth, task)
    pad = repmat('  ', 1, depth);
    cls = class(v);
    try sz = size(v); catch, sz = [0 0]; end
    szs = mat2str(sz);

    if isstruct(v) && isscalar(v)
        fn = fieldnames(v);
        fprintf(fid, '%s%s : struct  (%d fields)\n', pad, name, numel(fn));
        if depth >= 4, fprintf(fid, '%s  ...(max depth)\n', pad); return; end
        for i = 1:numel(fn)
            describe(fid, v.(fn{i}), fn{i}, depth+1, task);
        end
        return;
    end
    if isstruct(v)
        fprintf(fid, '%s%s : struct array %s\n', pad, name, szs);
        return;
    end
    if istable(v)
        vn = v.Properties.VariableNames;
        fprintf(fid, '%s%s : table  [%d rows x %d vars]\n', pad, name, height(v), width(v));
        for i = 1:numel(vn)
            col = v.(vn{i});
            extra = colSummary(col);
            fprintf(fid, '%s   %-16s %-10s %s\n', pad, vn{i}, class(col), extra);
        end
        % print first up-to-3 rows compactly
        nr = min(3, height(v));
        if nr > 0
            fprintf(fid, '%s   first %d row(s):\n', pad, nr);
            try
                txt = evalc('disp(head(v, nr))');
                txt = regexprep(txt, '\n', ['\n' pad '     ']);
                fprintf(fid, '%s     %s\n', pad, txt);
            catch
            end
        end
        return;
    end
    if iscell(v)
        fprintf(fid, '%s%s : cell %s\n', pad, name, szs);
        n = numel(v);
        showN = min(n, 12);
        items = {};
        for i = 1:showN
            items{end+1} = shortVal(v{i}); %#ok<AGROW>
        end
        if n > 0
            fprintf(fid, '%s     {%s%s}\n', pad, strjoin(items, ', '), tern(n>showN, ', ...', ''));
        end
        return;
    end
    if ischar(v)
        fprintf(fid, '%s%s : char  ''%s''\n', pad, name, oneLine(v));
        return;
    end
    if isstring(v)
        fprintf(fid, '%s%s : string %s  e.g. "%s"\n', pad, name, szs, oneLine(char(v(1))));
        return;
    end
    if isa(v, 'function_handle')
        fprintf(fid, '%s%s : function_handle  @%s\n', pad, name, func2str(v));
        return;
    end
    if islogical(v)
        fprintf(fid, '%s%s : logical %s  (%d true)\n', pad, name, szs, sum(v(:)));
        return;
    end
    if isnumeric(v)
        rng = '';
        try
            fv = double(v(:));
            fin = fv(isfinite(fv));
            if isempty(fin)
                rng = 'all NaN/Inf';
            else
                rng = sprintf('min=%.4g max=%.4g', min(fin), max(fin));
                nnan = sum(isnan(fv));
                if nnan>0, rng = sprintf('%s, %d NaN', rng, nnan); end
            end
        catch
        end
        samp = '';
        if numel(v) <= 8
            try samp = ['  vals=' mat2str(v, 5)]; catch, end
        end
        fprintf(fid, '%s%s : %s %s  %s%s\n', pad, name, cls, szs, rng, samp);
        return;
    end
    fprintf(fid, '%s%s : %s %s\n', pad, name, cls, szs);
end

function s = colSummary(col)
    s = '';
    try
        if isnumeric(col)
            fv = double(col(:)); fin = fv(isfinite(fv));
            if ~isempty(fin)
                s = sprintf('[%dx%d] min=%.4g max=%.4g', size(col,1), size(col,2), min(fin), max(fin));
                u = unique(fin);
                if numel(u) <= 8, s = [s '  uniq=' mat2str(u(:)',5)]; end
            else
                s = sprintf('[%dx%d]', size(col,1), size(col,2));
            end
        elseif iscell(col)
            ex = '';
            if ~isempty(col), ex = shortVal(col{1}); end
            s = sprintf('[%dx1 cell] e.g. %s', numel(col), ex);
        elseif iscategorical(col) || isstring(col)
            try
                u = unique(col); us = {};
                for i=1:min(8,numel(u)), us{end+1}=char(string(u(i))); end %#ok<AGROW>
                s = ['uniq={' strjoin(us, ',') '}'];
            catch, end
        end
    catch
    end
end

function s = shortVal(x)
    if ischar(x), s = ['''' oneLine(x) '''']; return; end
    if isstring(x), s = ['"' oneLine(char(x)) '"']; return; end
    if isnumeric(x) && isscalar(x), s = num2str(x); return; end
    if isnumeric(x), s = ['[' num2str(size(x,1)) 'x' num2str(size(x,2)) ' ' class(x) ']']; return; end
    s = class(x);
end

function s = oneLine(c)
    c = char(c); c = regexprep(c, '\s+', ' ');
    if numel(c) > 80, c = [c(1:77) '...']; end
    s = c;
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
