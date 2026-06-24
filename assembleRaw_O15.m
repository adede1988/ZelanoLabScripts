function raw = assembleRaw_O15(S)
%ASSEMBLERAW_O15  TASK-SPECIFIC raw load for O15 (raw_O15.mat + behavior CSV).
%   Does NOT run detect_ttls_O15 -- the dispatcher does that after the shared
%   common raw fields (type/paths) are set, because detection needs raw.paths.fig.

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
    if isstring(raw.labels), raw.labels = cellstr(raw.labels); end
end
