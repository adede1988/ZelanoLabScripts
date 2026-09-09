function build_session_index()
% BUILD_SESSION_INDEX  Master inventory of every sheet-listed session x task,
% cross-referenced against preproc finals on disk. Cheap: uses h5info only
% (no big-matrix loads). Writes out\tables\session_index.csv and prints a summary.
%
% Columns: task, sessID, cohort(Dupi/OBE/EEG), group(DupiS1/2/3/DupiSx/Control),
%          participant, sessNum, root, finalPath, onDisk, topVar, nChan, nSamp,
%          hasMoreThan1, moreThan1, hasBehDat, hasBmObj, hasBaseEmotion,
%          paramSource, missing(flag)

here = fileparts(mfilename('fullpath'));
proj = fileparts(here);                 % plotsForBruce_2
repo = fileparts(proj);                 % ZelanoLabScripts
addpath(repo);
outDir = fullfile(proj,'out','tables'); if ~exist(outDir,'dir'), mkdir(outDir); end

XLSX = 'R:\Neurology\Zelano_Lab\Lab_Common\Admin\Data\dataTracking.xlsx';

tasks   = {'breathingTask','cueTask','threshTask','O15'};
suffix  = containers.Map( ...
    {'breathingTask','cueTask','threshTask','O15'}, ...
    {'_breathingPreProc.mat','_cueTaskPreproc.mat','_PEA_threshold_preproc.mat','_O15preproc.mat'});
% fuzzy fallback globs
fuzzy   = containers.Map( ...
    {'breathingTask','cueTask','threshTask','O15'}, ...
    {'*breathing*.mat','*cueTask*.mat','*PEA_threshold*.mat','*O15*preproc*.mat'});

rows = {};
for ti = 1:numel(tasks)
    task = tasks{ti};
    fprintf('\n=== %s ===\n', task);
    cfg = applyParams(task,'main',XLSX);
    n = numel(cfg.sessionIDs);
    for r = 1:n
        id   = cfg.sessionIDs{r};
        root = cfg.root{r};
        ps   = ''; if ~isempty(cfg.paramSource), ps = cfg.paramSource{r}; end
        [cohort, group, participant, sessNum] = parse_cohort(id, cfg.datPrei(r));

        pdir = fullfile(root, id, 'preProc');
        fp   = fullfile(pdir, [id suffix(task)]);
        onDisk = exist(fp,'file')==2;
        if ~onDisk
            % fuzzy fallback
            dd = dir(fullfile(pdir, fuzzy(task)));
            dd = dd(~[dd.isdir]);
            if ~isempty(dd)
                % prefer a final ('preproc') over intermediate if both present
                nm = {dd.name};
                fp = fullfile(pdir, nm{1}); onDisk = true;
            end
        end

        topVar=""; nChan=NaN; nSamp=NaN; hasM1=false; m1=NaN;
        hasBeh=false; hasBm=false; hasBase=false;
        if onDisk
            try
                info = h5info(fp);
                % top-level group is the struct variable
                if ~isempty(info.Groups)
                    g = info.Groups(1);
                    topVar = string(strrep(g.Name,'/',''));
                    fnames = fieldNamesOf(g);
                    hasBeh  = any(strcmpi(fnames,'behDat'));
                    hasBm   = any(strcmpi(fnames,'bmObj'));
                    hasBase = any(strcmpi(fnames,'baseEmotion'));
                    hasM1   = any(strcmpi(fnames,'moreThan1'));
                    % data dims (MATLAB writes 2-D reversed: [nSamp nChan])
                    dg = findChild(g,'data');
                    if ~isempty(dg) && isfield(dg,'Dataspace') && ~isempty(dg.Dataspace)
                        sz = dg.Dataspace.Size;
                        if numel(sz)==2, nChan = min(sz); nSamp = max(sz); end
                    end
                    if hasM1
                        try, m1 = h5read(fp, [char(g.Name) '/moreThan1']); m1 = double(m1(1)); catch, end
                    end
                end
            catch e
                fprintf('  h5info FAIL %s : %s\n', id, e.message);
            end
        end

        missing = ~onDisk;
        rows(end+1,:) = {task, id, cohort, group, participant, sessNum, root, fp, ...
            double(onDisk), char(topVar), nChan, nSamp, double(hasM1), m1, ...
            double(hasBeh), double(hasBm), double(hasBase), ps, double(missing)}; %#ok<AGROW>
        if missing
            fprintf('  MISSING final: %s\n', id);
        end
    end
    fprintf('  %d sessions in sheet, %d missing on disk\n', n, sum(cell2mat(rows(strcmp(rows(:,1),task),19))));
end

T = cell2table(rows, 'VariableNames', {'task','sessID','cohort','group','participant', ...
    'sessNum','root','finalPath','onDisk','topVar','nChan','nSamp','hasMoreThan1', ...
    'moreThan1','hasBehDat','hasBmObj','hasBaseEmotion','paramSource','missing'});
writetable(T, fullfile(outDir,'session_index.csv'));
fprintf('\nWrote %s  (%d rows)\n', fullfile(outDir,'session_index.csv'), height(T));

% ---- quick summaries ----
fprintf('\n--- per task: onDisk / total ---\n');
for ti=1:numel(tasks)
    m = strcmp(T.task,tasks{ti});
    fprintf('  %-14s  %d / %d on disk\n', tasks{ti}, sum(T.onDisk(m)), sum(m));
end
fprintf('\n--- per cohort x task (onDisk) ---\n');
coh = {'Dupi','OBE','EEG'};
for c=1:3
  line = sprintf('  %-5s', coh{c});
  for ti=1:numel(tasks)
    m = strcmp(T.task,tasks{ti}) & strcmp(T.cohort,coh{c}) & T.onDisk==1;
    line = [line sprintf('  %-14s=%d', tasks{ti}, sum(m))]; %#ok<AGROW>
  end
  fprintf('%s\n', line);
end
end

% ================= helpers =================
function [cohort, group, participant, sessNum] = parse_cohort(id, datPrei)
    bits = strsplit(id,'_');
    participant = ''; sessNum = 1;
    if numel(bits)>=4, participant = bits{4}; end
    if numel(bits)>=5
        sn = str2double(bits{5});
        if ~isnan(sn), sessNum = sn; else, sessNum = 1; end
    else
        sessNum = 1;
    end
    % cohort from Type token (bits{2}) with datPrei as fallback
    typ = '';
    if numel(bits)>=2, typ = lower(bits{2}); end
    if contains(typ,'dupi'), cohort='Dupi';
    elseif contains(typ,'eeg'), cohort='EEG';
    elseif contains(typ,'obe'), cohort='OBE';
    else
        switch datPrei, case 1, cohort='Dupi'; case 3, cohort='EEG'; otherwise, cohort='OBE'; end
    end
    % TPB -> TB remap (documented quirk)
    if strcmpi(participant,'TPB'), participant='TB'; end
    if strcmp(cohort,'Dupi')
        if ismember(sessNum,[1 2 3]), group=sprintf('DupiS%d',sessNum); else, group='DupiSx'; end
    else
        group='Control';
    end
end

function fnames = fieldNamesOf(g)
    fnames = {};
    if isfield(g,'Datasets') && ~isempty(g.Datasets)
        fnames = [fnames, {g.Datasets.Name}];
    end
    if isfield(g,'Groups') && ~isempty(g.Groups)
        for k=1:numel(g.Groups)
            nm = g.Groups(k).Name; parts = strsplit(nm,'/'); fnames{end+1}=parts{end}; %#ok<AGROW>
        end
    end
end

function child = findChild(g, name)
    child = [];
    if isfield(g,'Datasets') && ~isempty(g.Datasets)
        i = find(strcmpi({g.Datasets.Name}, name),1);
        if ~isempty(i), child = g.Datasets(i); return; end
    end
    if isfield(g,'Groups') && ~isempty(g.Groups)
        for k=1:numel(g.Groups)
            parts = strsplit(g.Groups(k).Name,'/');
            if strcmpi(parts{end}, name), child = g.Groups(k); return; end
        end
    end
end
