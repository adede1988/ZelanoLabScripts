function o15_explore()
% Discovery + structure exploration for the O15 preprocessed finals.
% Lists which O15 finals exist on disk (+ mod time, flags today's), then loads
% the most-recently-modified one and dumps its structure (TTL, behDat, labels).

repo = fileparts(fileparts(mfilename('fullpath')));   % ...\ZelanoLabScripts
addpath(repo);

cfg = applyParams('O15','main');
ids   = cfg.sessionIDs;
roots = cfg.root;
n = numel(ids);
fprintf('applyParams O15 main: %d sessions listed\n', n);

todayNum = datenum(2026,6,24);
existing = {};  % {id, path, datenum, isToday}
for i = 1:n
    id = ids{i};
    f = fullfile(roots{i}, id, 'preProc', [id '_O15preproc.mat']);
    if ~isfile(f)
        d = dir(fullfile(roots{i}, id, 'preProc', [id '*O15*preproc*.mat']));
        if ~isempty(d), f = fullfile(d(1).folder, d(1).name); end
    end
    if isfile(f)
        dn = dir(f); dn = dn(1).datenum;
        isToday = dn >= todayNum;
        existing(end+1,:) = {id, f, dn, isToday}; %#ok<AGROW>
    end
end
fprintf('\n==== O15 finals on disk: %d ====\n', size(existing,1));
for i = 1:size(existing,1)
    fprintf('  %s  | %s | %s\n', existing{i,1}, datestr(existing{i,3},'yyyy-mm-dd HH:MM'), ...
        ternary(existing{i,4},'TODAY',''));
end
if ~isempty(existing)
    nToday = sum([existing{:,4}]);
    fprintf('  -> modified TODAY (6/24): %d\n', nToday);
end

if isempty(existing), fprintf('NO O15 finals found.\n'); return; end

% pick the most-recently-modified final to explore
[~, mi] = max([existing{:,3}]);
fp = existing{mi,2}; id = existing{mi,1};
fprintf('\n==== EXPLORING %s ====\n%s\n', id, fp);

S = load(fp); fn = fieldnames(S); od = S.(fn{1});
fprintf('top-level var: %s\n', fn{1});
fprintf('outDat fields:\n'); disp(fieldnames(od)');

fprintf('\n-- basics --\n');
fprintf('task=%s type=%s fs=%g sessID=%s\n', getf(od,'task'), getf(od,'type'), od.fs, getf(od,'sessID'));
fprintf('data size: [%d x %d]\n', size(od.data,1), size(od.data,2));
fprintf('moreThan1=%s\n', num2str(getfo(od,'moreThan1')));

fprintf('\n-- labels --\n');
labs = cellfun(@(x) char(string(x)), od.labels, 'uni', 0);
for i = 1:numel(labs), fprintf('  %2d: %s\n', i, labs{i}); end
isMac = cellfun(@(x) contains(x,'macBP'), labs);
fprintf('macBP channels: %s\n', strjoin(labs(isMac), ', '));
isR = cellfun(@(x) contains(x,'rsp'), labs);
fprintf('rsp channels: %s | rspIDX=%g rspFlip=%g\n', strjoin(labs(isR), ', '), getfo(od,'rspIDX'), getfo(od,'rspFlip'));

fprintf('\n-- TTL --\n');
if isfield(od,'TTL')
    T = od.TTL;
    fprintf('class: %s\n', class(T));
    if istable(T)
        fprintf('size: [%d x %d]\n', height(T), width(T));
        fprintf('vars: %s\n', strjoin(T.Properties.VariableNames, ', '));
        fprintf('head:\n'); disp(head(T, 5));
    else
        fprintf('numel=%d, first 10: %s\n', numel(T), mat2str(T(1:min(10,numel(T)))));
    end
end

fprintf('\n-- behDat --\n');
if isfield(od,'behDat')
    bd = od.behDat;
    fprintf('class: %s\n', class(bd));
    if istable(bd)
        fprintf('size: [%d x %d]\n', height(bd), width(bd));
        fprintf('vars: %s\n', strjoin(bd.Properties.VariableNames, ', '));
        fprintf('head(8):\n'); disp(head(bd, 8));
        if ismember('sniffLabel', bd.Properties.VariableNames)
            sl = string(bd.sniffLabel);
            u = unique(sl);
            fprintf('sniffLabel categories + counts:\n');
            for k = 1:numel(u), fprintf('   %-10s : %d\n', u(k), sum(sl==u(k))); end
        end
        for vv = {'wiTriali','n','sniffType'}
            if ismember(vv{1}, bd.Properties.VariableNames)
                x = bd.(vv{1});
                if isnumeric(x), fprintf('%s range: [%g .. %g], unique=%s\n', vv{1}, min(x), max(x), mat2str(unique(x(:))')); end
            end
        end
        for vv = {'sniffOnset','finalOnset','adjust','TTLoffSet'}
            if ismember(vv{1}, bd.Properties.VariableNames)
                x = bd.(vv{1});
                fprintf('%s: min=%g max=%g (first 6: %s)\n', vv{1}, min(x), max(x), mat2str(x(1:min(6,numel(x)))'));
            end
        end
        if ismember('sniffLabel',bd.Properties.VariableNames) && ismember('finalOnset',bd.Properties.VariableNames)
            isStart = string(bd.sniffLabel)=="start";
            fo = bd.finalOnset;
            fprintf('\nfirst few START sniff finalOnsets: %s\n', mat2str(fo(find(isStart,5))'));
            fprintf('first START finalOnset=%g (%.1f s into recording); recording length=%.1f s\n', ...
                min(fo(isStart)), min(fo(isStart))/od.fs, size(od.data,2)/od.fs);
        end
    end
end

for vv = {'CSClist','OGdataDir','loadFile','preProcScript'}
    if isfield(od, vv{1})
        v = od.(vv{1});
        if ischar(v) || isstring(v), fprintf('%s = %s\n', vv{1}, char(string(v)));
        else, fprintf('%s present (class %s)\n', vv{1}, class(v)); end
    end
end
fprintf('\n==== done ====\n');
end

function s = getf(od,f), if isfield(od,f), s=char(string(od.(f))); else, s='<none>'; end, end
function v = getfo(od,f), if isfield(od,f) && ~isempty(od.(f)), v=od.(f); else, v=NaN; end, end
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
