function T = o15_session_table(verbose)
% O15_SESSION_TABLE  Cheap manifest of O15-task sessions (no multi-GB loads).
%   T = o15_session_table()  returns one row per O15 session listed by
%   applyParams('O15','main'), with identity/group info, on-disk status, cheap
%   file metadata (channel/sample counts via h5info) and a `fresh` flag.
%
%   `fresh` = the final .mat was modified on/after 2026-06-24 (the date O15
%   preprocessing began on this run). Per request we ONLY analyse fresh finals;
%   older finals on disk (a previous preprocessing pass) are listed but excluded
%   by the drivers (which filter T(T.onDisk & T.fresh,:)).
%
%   Group rule (matches the cue analysis): DupiS1/S2/S3 = type Dupi & trailing
%   _1/_2/_3; Control = all non-Dupi; Dupi with sessNum>3 -> 'DupiSx' (reported
%   but excluded from the 4-group means).

    if nargin < 1, verbose = true; end
    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo);

    freshCutoff = datenum(2026,6,24);    % O15 preprocessing (re)started this date

    cfg = applyParams('O15','main');
    ids   = cfg.sessionIDs;
    roots = cfg.root;
    n = numel(ids);

    sessID = strings(n,1); subID = strings(n,1); typ = strings(n,1);
    sessNum = zeros(n,1);  grp = strings(n,1);
    fpath = strings(n,1);  exists = false(n,1); fresh = false(n,1); modWhen = strings(n,1);
    nChan = nan(n,1); nSamp = nan(n,1);

    for i = 1:n
        id = ids{i};
        sessID(i) = string(id);
        f = fullfile(roots{i}, id, 'preProc', [id '_O15preproc.mat']);
        if ~isfile(f)
            d = dir(fullfile(roots{i}, id, 'preProc', [id '*O15*preproc*.mat']));
            if ~isempty(d), f = fullfile(d(1).folder, d(1).name); end
        end
        fpath(i) = string(f);
        exists(i) = isfile(f);
        if exists(i)
            di = dir(f);
            fresh(i) = di(1).datenum >= freshCutoff;
            modWhen(i) = string(datestr(di(1).datenum,'yyyy-mm-dd HH:MM'));
        end

        bits = strsplit(id,'_');
        if numel(bits) >= 4, subID(i) = string(bits{4}); end
        if numel(bits) >= 5, sessNum(i) = str2double(bits{5}); else, sessNum(i) = 1; end
        try
            P = applyParams('O15', id);
            typ(i) = string(P.type);
        catch
            typ(i) = "";
        end
        if strcmpi(typ(i),'Dupi')
            switch sessNum(i)
                case 1, grp(i) = "DupiS1";
                case 2, grp(i) = "DupiS2";
                case 3, grp(i) = "DupiS3";
                otherwise, grp(i) = "DupiSx";
            end
        else
            grp(i) = "Control";
        end

        if exists(i)
            try
                info = h5info(f, '/outDat/data');
                sz = info.Dataspace.Size;   % MATLAB writes 2-D as [nSamp nChan]
                nSamp(i) = sz(1); nChan(i) = sz(2);
            catch
                for vn = ["out","chanDat"]
                    try
                        info = h5info(f, ['/' char(vn) '/data']);
                        sz = info.Dataspace.Size; nSamp(i)=sz(1); nChan(i)=sz(2); break;
                    catch
                    end
                end
            end
        end
    end

    T = table(sessID, subID, sessNum, typ, grp, exists, fresh, modWhen, nChan, nSamp, fpath, ...
        'VariableNames', {'sessID','subID','sessNum','type','group','onDisk','fresh','modWhen','nChan','nSamp','path'});

    if verbose
        Td = T(T.onDisk & T.fresh,:);
        fprintf('\nO15 sessions listed: %d | on disk: %d | fresh (>=6/24): %d\n', ...
            height(T), sum(T.onDisk), height(Td));
        fprintf('Group counts (fresh on disk):\n');
        g = categorical(Td.group);
        cats = categories(g);
        for k = 1:numel(cats)
            fprintf('   %-8s : %d\n', cats{k}, sum(g==cats{k}));
        end
        fprintf('\n%s\n', evalc('disp(T)'));
    end
end
