function outDat = assemble_outDat_breathingTask(raw, S, P)
% Mirror O15's assembler so the rest of the pipeline looks identical.

outDat = struct();
outDat.behDat   = raw.beh;
outDat.labels   = raw.labels;
outDat.fs       = raw.fs_raw;
outDat.data     = raw.data;
outDat.sessID   = S.id;
outDat.task     = P.task;
outDat.figs     = fullfile(S.fig, outDat.task);
if ~exist(outDat.figs,'dir'), mkdir(outDat.figs); end
outDat.type     = P.type;

% carry through any original TTL container for the cue task
if isfield(raw,'TTL'), outDat.TTL = raw.TTL; end

% optional placeholders (keeps parity with O15 fields)
outDat.rspIDX   = P.rspIDX;
outDat.rspFlip  = P.rspFlip;

end
