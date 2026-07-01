function T = chirp_session_table(C, freshOnly, verbose)
% CHIRP_SESSION_TABLE  Discover cue-task preprocessed finals on E: (no R: needed).
%   T = chirp_session_table(C, freshOnly, verbose)
%
%   Scans C.dataRoot (default E:\Lab_Common) for <id>\preProc\<id>_cueTaskPreproc.mat,
%   derives identity/group like cue_session_table, and flags each file's freshness against
%   C.freshCutoff (D5: the rerun's corrected bestMac/noise only exist in finals modified
%   after 4pm 2026-06-29).  Cheap: file metadata only, no .mat loads.
%
%   freshOnly (default true) -> return only finals with mtime > C.freshCutoff.
%
%   Columns: sessID subID sessNum type group path mtime isFresh nMacBP_hint
%   (nMacBP_hint via h5info, best-effort; NaN if unreadable).

    if nargin < 1 || isempty(C), C = chirp_config(); end
    if nargin < 2 || isempty(freshOnly), freshOnly = true; end
    if nargin < 3 || isempty(verbose), verbose = true; end

    root = C.dataRoot;
    d = dir(fullfile(root, '**', 'preProc', '*_cueTaskPreproc.mat'));
    % case/naming fallback: also catch *_cueTaskPreProc.mat (older intermediate name)
    d2 = dir(fullfile(root, '**', 'preProc', '*cueTask*reproc.mat'));
    d = dedupeFiles([d; d2]);

    n = numel(d);
    sessID=strings(n,1); subID=strings(n,1); typ=strings(n,1); grp=strings(n,1);
    sessNum=zeros(n,1); fpath=strings(n,1); mtime=NaT(n,1); isFresh=false(n,1);
    nMac=nan(n,1);

    cutoff = C.freshCutoff;
    for i = 1:n
        f = fullfile(d(i).folder, d(i).name);
        fpath(i) = string(f);
        mtime(i) = datetime(d(i).datenum, 'ConvertFrom','datenum');
        isFresh(i) = mtime(i) > cutoff;

        % sessID = folder name two levels up (...\<id>\preProc\<file>)
        [pp,~] = fileparts(d(i).folder);     % ...\<id>
        [~,id] = fileparts(pp);
        sessID(i) = string(id);

        bits = strsplit(id, '_');
        if numel(bits) >= 4, subID(i) = string(bits{4}); end
        if numel(bits) >= 5, sessNum(i) = str2double(bits{5}); else, sessNum(i) = 1; end
        if isnan(sessNum(i)), sessNum(i) = 1; end

        % type from the path root (...\Lab_Common\Dupi\...  vs  ...\OBEControl\...)
        fl = lower(f);
        if contains(fl, [filesep 'dupi' filesep]),        typ(i) = "Dupi";
        elseif contains(fl, [filesep 'obecontrol' filesep]), typ(i) = "OBE";
        elseif numel(bits) >= 2,                            typ(i) = string(bits{2});
        else,                                               typ(i) = ""; end

        if strcmpi(typ(i), 'Dupi')
            switch sessNum(i)
                case 1, grp(i) = "DupiS1"; case 2, grp(i) = "DupiS2";
                case 3, grp(i) = "DupiS3"; otherwise, grp(i) = "DupiSx";
            end
        else
            grp(i) = "Control";
        end

        % cheap macBP-count hint via h5info (HDF5 -v7.3); best-effort
        try
            info = h5info(f, '/outDat/labels');  %#ok<NASGU>  % presence check only
            nMac(i) = NaN; % label strings not trivially countable from h5info; leave to loader
        catch
        end
    end

    T = table(sessID, subID, sessNum, typ, grp, fpath, mtime, isFresh, nMac, ...
        'VariableNames', {'sessID','subID','sessNum','type','group','path','mtime','isFresh','nMacBP_hint'});
    T = sortrows(T, {'isFresh','mtime'}, {'descend','ascend'});

    if freshOnly, T = T(T.isFresh, :); end

    if verbose
        fprintf('chirp_session_table: %d cue finals under %s\n', n, root);
        fprintf('  fresh (> %s): %d   stale: %d\n', char(cutoff), sum(isFresh), sum(~isFresh));
        if freshOnly, fprintf('  returning FRESH only: %d rows\n', height(T)); end
        if ~isempty(T)
            g = categorical(T.group); cats = categories(g);
            for k = 1:numel(cats), fprintf('    %-8s : %d\n', cats{k}, sum(g==cats{k})); end
        end
    end
end

function d = dedupeFiles(d)
    full = arrayfun(@(x) lower(fullfile(x.folder, x.name)), d, 'uni', 0);
    [~, ia] = unique(full, 'stable');
    d = d(ia);
end
