function raw = assembleRaw_breathingTask(S)
%ASSEMBLERAW_BREATHINGTASK  TASK-SPECIFIC raw load for breathingTask.
%   Loads <id>_breathingPreProc.mat (shared loadIntermediateRaw) and applies the
%   breathing-specific TTL handling: a 5-min-window fallback when the file has no
%   TTL, then the /4 sample-rate rescale. Called by the assemble_outDat_all
%   dispatcher; shared assembly is in assembleOutDat.

    matPath = fullfile(S.root, S.id, 'preProc', [S.id '_breathingPreProc.mat']);
    [raw, od] = loadIntermediateRaw(matPath);

    if ~isfield(od, 'TTL')
        % 5-min window fallback (assumes contiguous concatenated blocks)
        TTLv = 0:600000:size(od.data, 2);
        TTLv(1)   = 1;
        TTLv(end) = [];
        od.TTL = TTLv;
    end
    raw.TTL = round(od.TTL ./ 4);
end
