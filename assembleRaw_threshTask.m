function raw = assembleRaw_threshTask(S)
%ASSEMBLERAW_THRESHTASK  TASK-SPECIFIC raw load for threshTask.
%   Loads <id>_PEA_threshold_preproc.mat (shared loadIntermediateRaw) and passes
%   through its TTL table. Called by the assemble_outDat_all dispatcher.

    matPath = fullfile(S.root, S.id, 'preProc', [S.id '_PEA_threshold_preproc.mat']);
    [raw, od] = loadIntermediateRaw(matPath);
    if isfield(od, 'TTL'), raw.TTL = od.TTL; end
end
