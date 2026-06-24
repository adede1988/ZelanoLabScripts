function outDat = assembleOutDat(raw, S, P, figDir)
%ASSEMBLEOUTDAT  SHARED assembly of the common outDat struct from raw.
%
%   outDat = assembleOutDat(raw, S, P, figDir)
%
%   No task-specific branches. The only conditional is the data-driven TTL
%   passthrough: raw carries a .TTL field for breathing/cue/thresh but not for
%   O15 (whose TTL is returned separately by detect_ttls_O15 and attached by the
%   O15 main), so `isfield(raw,'TTL')` is the correct, task-agnostic gate.
%   O15-specific outDat fields are added by assembleOutDat_O15extras.

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

    if isfield(raw, 'TTL')
        outDat.TTL = raw.TTL;
    end
end
