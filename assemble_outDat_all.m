function [outDat, raw, TTL] = assemble_outDat_all(S, P)
% assemble_outDat_all  BACKWARD-COMPATIBLE wrapper (one-call assemble API).
%
%   [outDat, raw, TTL] = assemble_outDat_all(S, P)
%
%   The deliverable *_main.m scripts do NOT use this -- they call the
%   task-specific loader + shared assembler directly, so the task-specific
%   pieces are visible right in each main. This wrapper preserves the old
%   one-call API for the _dev/run_* harnesses and any other callers, and simply
%   composes the same functions:
%     task-specific : assembleRaw_breathingTask / _cueTask / _threshTask / _O15,
%                     assembleOutDat_O15extras
%     shared        : assembleOutDat

    switch char(P.task)
        case 'breathingTask', raw = assembleRaw_breathingTask(S); TTL = [];
        case 'cueTask',       raw = assembleRaw_cueTask(S);       TTL = [];
        case 'threshTask',    raw = assembleRaw_threshTask(S);    TTL = [];
        case 'O15',           [raw, TTL] = assembleRaw_O15(S, P);
        otherwise
            error('assemble_outDat_all:badTask', 'Unsupported task "%s".', char(P.task));
    end

    outDat = assembleOutDat(raw, S, P);

    if strcmp(char(P.task), 'O15')
        outDat = assembleOutDat_O15extras(outDat, S, raw);
    end
end
