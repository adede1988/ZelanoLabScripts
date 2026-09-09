function T = thresh_session_table(verbose)
% THRESH_SESSION_TABLE  Cheap manifest of thresh-task sessions (no multi-GB loads).
%   T = thresh_session_table()  returns a table with one row per thresh session
%   that has a final preproc .mat on disk, with identity/group info and cheap file
%   metadata (channel count, sample count) read via h5info.
%
%   Final filename: <id>_PEA_threshold_preproc.mat  (top-level var 'outDat').
%
%   Group rule (4 groups): DupiS1/DupiS2/DupiS3 = type Dupi & trailing _1/_2/_3;
%   Control = all non-Dupi. Dupi with sessNum>3 -> group 'DupiSx' (excluded from
%   the 4-group means).

    if nargin < 1, verbose = true; end
    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(repo);

    cfg = applyParams('threshTask','main');
    ids   = cfg.sessionIDs;
    roots = cfg.root;
    n = numel(ids);

    sessID = strings(n,1); subID = strings(n,1); typ = strings(n,1);
    sessNum = zeros(n,1);  grp = strings(n,1);
    fpath = strings(n,1);  exists = false(n,1);
    nChan = nan(n,1); nSamp = nan(n,1);

    for i = 1:n
        id = ids{i};
        sessID(i) = string(id);
        f = fullfile(roots{i}, id, 'preProc', [id '_PEA_threshold_preproc.mat']);
        if ~isfile(f)
            % case-insensitive / fuzzy fallback
            d = dir(fullfile(roots{i}, id, 'preProc', [id '*PEA_threshold*.mat']));
            if ~isempty(d), f = fullfile(d(1).folder, d(1).name); end
        end
        fpath(i) = string(f);
        exists(i) = isfile(f);

        % identity from sessID + sheet
        bits = strsplit(id,'_');
        if numel(bits) >= 4, subID(i) = string(bits{4}); end
        if numel(bits) >= 5, sessNum(i) = str2double(bits{5}); else, sessNum(i) = 1; end
        try
            P = applyParams('threshTask', id);
            typ(i) = string(P.type);
        catch
            typ(i) = "";
        end
        % group
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

        % cheap file metadata via h5info (no data load)
        if exists(i)
            try
                info = h5info(f, '/outDat/data');
                sz = info.Dataspace.Size;   % HDF5 dim order; MATLAB [nChan nSamp]
                % MATLAB writes 2-D as [nSamp nChan] in the file -> reverse
                nSamp(i) = sz(1); nChan(i) = sz(2);
            catch
                % some files store var as 'out'/'chanDat'
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

    T = table(sessID, subID, sessNum, typ, grp, exists, nChan, nSamp, fpath, ...
        'VariableNames', {'sessID','subID','sessNum','type','group','onDisk','nChan','nSamp','path'});

    if verbose
        Td = T(T.onDisk,:);
        fprintf('\nThresh sessions listed: %d | on disk: %d\n', height(T), height(Td));
        fprintf('Group counts (on disk):\n');
        g = categorical(Td.group);
        cats = categories(g);
        for k = 1:numel(cats)
            fprintf('   %-8s : %d\n', cats{k}, sum(g==cats{k}));
        end
        fprintf('\n%s\n', evalc('disp(T)'));
    end
end
