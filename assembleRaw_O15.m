function [raw, TTL] = assembleRaw_O15(S, P)
%ASSEMBLERAW_O15  TASK-SPECIFIC raw load for O15.
%   Loads raw_O15.mat + the behavior CSV, then parses the photodiode into the
%   TTL table via detect_ttls_O15 (which saves TTLs.jpg into the figure dir).
%   Returns the raw struct and the O15 TTL table. The O15-only outDat fields are
%   added afterwards by assembleOutDat_O15extras.

    matPath = fullfile(S.root, S.id, 'raw', 'raw_O15', 'raw_O15.mat');
    if ~exist(matPath, 'file')
        error('assembleRaw_O15:MissingMat', 'Raw MAT not found: %s', matPath);
    end
    dat = load(matPath);
    if isfield(dat, 'curDat')
        dat = dat.curDat;
    else
        error('assembleRaw_O15:BadMat', 'Expected curDat in %s', matPath);
    end

    behPath = fullfile(S.root, S.id, 'Behavioral_data', 'O15', ...
                       sprintf('O15_responses_%s.csv', S.id));
    if ~exist(behPath, 'file')
        error('assembleRaw_O15:MissingBehavior', 'Behavior CSV not found: %s', behPath);
    end

    raw.sessID = char(S.id);
    raw.fs_raw = dat.rawData.fsample;
    raw.data   = dat.rawData.trial{1};
    raw.labels = dat.outLabs;
    if isfield(dat, 'ncslabels'), raw.ncslabels = dat.ncslabels; end
    raw.beh    = readtable(behPath);
    raw.paths.mat = matPath;
    raw.paths.beh = behPath;
    raw.paths.fig = resolveFigDir(S);   % detect_ttls_O15 saves TTLs.jpg here
    if isstring(raw.labels), raw.labels = cellstr(raw.labels); end

    % O15-specific: photodiode -> TTL table
    [TTL, raw] = detect_ttls_O15(raw, P);
end
