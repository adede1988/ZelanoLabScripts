function [outDat, raw, TTL] = assemble_outDat_all(S, P)
% assemble_outDat_all  Combined raw-load + outDat assembly for all four tasks.
%
%   [outDat, raw, TTL] = assemble_outDat_all(S, P)
%
%   Replaces (getSessionParams_<task> raw-load half) + (detect_ttls_O15 for O15)
%   + (assemble_outDat_<task>), branching on P.task.
%
%   Inputs
%     S : struct with .id, .root, .fig  (and .figPath for O15).
%     P : parameter struct from applyParams (P.task selects the branch).
%
%   Outputs
%     outDat : assembled struct fed into the shared downstream pipeline.
%     raw    : raw view (used after assemble by build_behavior_table_* in
%              cue/thresh/O15).
%     TTL    : O15 TTL table from detect_ttls_O15; [] for the other three tasks.

    task = char(P.task);

    % --- figure dir ---
    if isfield(S, 'fig') && ~isempty(S.fig)
        figDir = S.fig;
    elseif isfield(S, 'figPath') && ~isempty(S.figPath)
        figDir = fullfile(S.figPath, S.id);
    else
        error('assemble_outDat_all:noFigDir', 'S needs .fig or .figPath.');
    end
    if ~isfolder(figDir), mkdir(figDir); end

    raw = struct();
    TTL = [];

    % ================= TASK-SPECIFIC: raw load per task =================
    % Each case knows where a task's raw data lives and how to read it (and,
    % for O15, runs detect_ttls_O15). To add a new task, add a new case here.
    % Everything AFTER the switch (common raw fields + outDat assembly) is SHARED.
    switch task
        % =================================================================
        case 'breathingTask'
            matPath = fullfile(S.root, S.id, 'preProc', [S.id '_breathingPreProc.mat']);
            if ~exist(matPath, 'file')
                error('assemble_outDat_all:MissingMat', ...
                    'Expected breathing MAT at %s.', matPath);
            end
            tmp = load(matPath);
            if ~isfield(tmp, 'outDat') && ~isfield(tmp, 'out') && ~isfield(tmp, 'chanDat')
                error('assemble_outDat_all:NoOutDat', 'MAT must contain outDat/out/chanDat.');
            end
            try
                try
                    od = tmp.outDat;
                catch
                    od = tmp.chanDat;
                end
            catch
                od = tmp.out;
            end

            raw.sessID = char(od.sessID);
            raw.fs_raw = od.fs;
            raw.data   = od.data;
            raw.labels = od.labels;
            raw.beh    = od.behDat;

            if ~isfield(od, 'TTL')
                % 5-min window fallback (assumes contiguous concatenated blocks)
                TTLv = 0:600000:size(od.data, 2);
                TTLv(1)   = 1;
                TTLv(end) = [];
                od.TTL = TTLv;
            end
            raw.TTL = round(od.TTL ./ 4);

        % =================================================================
        case 'cueTask'
            matPath = fullfile(S.root, S.id, 'preProc', [S.id '_cueTaskPreProc.mat']);
            if ~exist(matPath, 'file')
                error('assemble_outDat_all:MissingMat', ...
                    'Expected cue MAT at %s.', matPath);
            end
            tmp = load(matPath);
            if ~isfield(tmp, 'outDat')
                error('assemble_outDat_all:NoOutDat', 'MAT must contain outDat.');
            end
            od = tmp.outDat;

            raw.sessID = char(od.sessID);
            raw.fs_raw = od.fs;
            raw.data   = od.data;
            raw.labels = od.labels;
            raw.beh    = od.behDat;
            if isfield(od, 'TTL'), raw.TTL = od.TTL; end

        % =================================================================
        case 'threshTask'
            matPath = fullfile(S.root, S.id, 'preProc', [S.id '_PEA_threshold_preproc.mat']);
            if ~exist(matPath, 'file')
                error('assemble_outDat_all:MissingMat', ...
                    'Expected thresh MAT at %s.', matPath);
            end
            tmp = load(matPath);
            if ~isfield(tmp, 'outDat')
                error('assemble_outDat_all:NoOutDat', 'MAT must contain outDat.');
            end
            od = tmp.outDat;

            raw.sessID = char(od.sessID);
            raw.fs_raw = od.fs;
            raw.data   = od.data;
            raw.labels = od.labels;
            raw.beh    = od.behDat;
            if isfield(od, 'TTL'), raw.TTL = od.TTL; end

        % =================================================================
        case 'O15'
            matPath = fullfile(S.root, S.id, 'raw', 'raw_O15', 'raw_O15.mat');
            if ~exist(matPath, 'file')
                error('assemble_outDat_all:MissingMat', 'Raw MAT not found: %s', matPath);
            end
            dat = load(matPath);
            if isfield(dat, 'curDat')
                dat = dat.curDat;
            else
                error('assemble_outDat_all:BadMat', 'Expected curDat in %s', matPath);
            end

            behPath = fullfile(S.root, S.id, 'Behavioral_data', 'O15', ...
                               sprintf('O15_responses_%s.csv', S.id));
            if ~exist(behPath, 'file')
                error('assemble_outDat_all:MissingBehavior', 'Behavior CSV not found: %s', behPath);
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

        otherwise
            error('assemble_outDat_all:badTask', 'Unsupported task "%s".', task);
    end

    % --- common raw fields ---
    raw.type = P.type;
    raw.paths.root = S.root;
    raw.paths.fig  = figDir;

    % --- O15: TTL detection (needs raw.paths.fig for the TTLs.jpg) ---
    if strcmp(task, 'O15')
        [TTL, raw] = detect_ttls_O15(raw, P);
    end

    % --- assemble outDat ---
    outDat = struct();
    outDat.behDat = raw.beh;
    outDat.labels = raw.labels;
    outDat.fs     = raw.fs_raw;
    outDat.data   = raw.data;
    outDat.sessID = S.id;
    outDat.task   = P.task;
    outDat.type   = P.type;
    outDat.figs   = fullfile(figDir, char(P.task));
    if ~exist(outDat.figs, 'dir'), mkdir(outDat.figs); end
    outDat.rspIDX  = P.rspIDX;
    outDat.rspFlip = P.rspFlip;

    if strcmp(task, 'O15')
        outDat.CSClist = raw.ncslabels;
        outDat.OGdataDir = fullfile(S.root, S.id);
        tmp = dir(fullfile(S.root, S.id));
        tmp = tmp(cellfun(@(x) contains(x, '.m'), {tmp.name}));
        tmp = tmp(cellfun(@(x) contains(x, 'LoadData'), {tmp.name}));
        if size(tmp, 1) == 1
            outDat.loadFile = tmp.name;
        else
            error('assemble_outDat_all:loadFile', 'load file not identified uniquely');
        end
        outDat.preProcScript = 'O15PreProc.m';
        % NB: do NOT set outDat.TTL here for O15 -- the main script does that.
    else
        if isfield(raw, 'TTL')
            outDat.TTL = raw.TTL;
        end
    end
end
