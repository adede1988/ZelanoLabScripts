function [raw, od] = loadIntermediateRaw(matPath)
%LOADINTERMEDIATERAW  SHARED loader for a <task>PreProc.mat intermediate.
%
%   [raw, od] = loadIntermediateRaw(matPath)
%
%   Loads the intermediate, extracts the stored struct (saved as outDat, or
%   chanDat via parSave, or out in older files), and sets the common raw fields.
%   Returns the loaded struct `od` too, so the caller can apply its own
%   task-specific TTL handling. No task-specific branches.

    if ~exist(matPath, 'file')
        error('loadIntermediateRaw:MissingMat', 'Expected intermediate MAT at %s.', matPath);
    end
    tmp = load(matPath);
    if ~isfield(tmp, 'outDat') && ~isfield(tmp, 'out') && ~isfield(tmp, 'chanDat')
        error('loadIntermediateRaw:NoOutDat', 'MAT must contain outDat/out/chanDat: %s', matPath);
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
end
