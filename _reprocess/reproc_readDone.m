function ids = reproc_readDone(p)
% Read the set of already-reprocessed session IDs from a done-marker file.
    ids = {};
    if exist(p, 'file') == 2
        v = strtrim(string(splitlines(fileread(p))));
        ids = cellstr(v(strlength(v) > 0));
    end
end
