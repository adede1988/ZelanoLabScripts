function [outDat, raw, TTL] = assemble_outDat_all(S, P)
% assemble_outDat_all  DISPATCHER: route to the task-specific raw loader, set the
%   shared common raw fields, then build outDat with the shared assembler.
%
%   [outDat, raw, TTL] = assemble_outDat_all(S, P)
%
%   This is the single place in the pipeline that branches on task. All actual
%   work lives in functions that are EITHER completely task-specific OR completely
%   shared:
%     task-specific : assembleRaw_breathingTask / _cueTask / _threshTask / _O15,
%                     detect_ttls_O15, assembleOutDat_O15extras
%     shared        : assembleOutDat, resolveFigDir (below)
%
%   Inputs
%     S : struct with .id, .root, .fig  (and .figPath for O15).
%     P : parameter struct from applyParams (P.task selects the loader).
%   Outputs
%     outDat : assembled struct fed into the shared downstream pipeline.
%     raw    : raw view (used after assemble by build_behavior_table_* in
%              cue/thresh/O15).
%     TTL    : O15 TTL table from detect_ttls_O15; [] for the other three tasks.

    task   = char(P.task);
    figDir = resolveFigDir(S);
    TTL    = [];

    % --- route to the task-specific raw loader ---
    switch task
        case 'breathingTask', raw = assembleRaw_breathingTask(S);
        case 'cueTask',       raw = assembleRaw_cueTask(S);
        case 'threshTask',    raw = assembleRaw_threshTask(S);
        case 'O15',           raw = assembleRaw_O15(S);
        otherwise
            error('assemble_outDat_all:badTask', 'Unsupported task "%s".', task);
    end

    % --- shared: common raw fields ---
    raw.type       = P.type;
    raw.paths.root = S.root;
    raw.paths.fig  = figDir;

    % --- O15 only: photodiode -> TTL table (needs raw.paths.fig for TTLs.jpg) ---
    if strcmp(task, 'O15')
        [TTL, raw] = detect_ttls_O15(raw, P);
    end

    % --- shared: assemble the common outDat ---
    outDat = assembleOutDat(raw, S, P, figDir);

    % --- O15 only: extra outDat fields ---
    if strcmp(task, 'O15')
        outDat = assembleOutDat_O15extras(outDat, S, raw);
    end
end

% ============================ shared helper ============================

function figDir = resolveFigDir(S)
    if isfield(S, 'fig') && ~isempty(S.fig)
        figDir = S.fig;
    elseif isfield(S, 'figPath') && ~isempty(S.figPath)
        figDir = fullfile(S.figPath, S.id);
    else
        error('assemble_outDat_all:noFigDir', 'S needs .fig or .figPath.');
    end
    if ~isfolder(figDir), mkdir(figDir); end
end
