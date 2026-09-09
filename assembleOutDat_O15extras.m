function outDat = assembleOutDat_O15extras(outDat, S, raw)
%ASSEMBLEOUTDAT_O15EXTRAS  TASK-SPECIFIC O15 outDat fields.
%   Adds the CSC list, the original-data directory, the session's LoadData
%   script (must be unique), and the preproc-script name. O15 only.

    outDat.CSClist   = raw.ncslabels;
    outDat.OGdataDir = fullfile(S.root, S.id);

    tmp = dir(fullfile(S.root, S.id));
    tmp = tmp(cellfun(@(x) contains(x, '.m'), {tmp.name}));
    tmp = tmp(cellfun(@(x) contains(x, 'LoadData'), {tmp.name}));
    if size(tmp, 1) == 1
        outDat.loadFile = tmp.name;
    else
        tmp = tmp(cellfun(@(x) contains(x, 'AD'), {tmp.name}));
         if size(tmp, 1) == 1
            outDat.loadFile = tmp.name;
        else
            error('assembleOutDat_O15extras:loadFile', 'load file not identified uniquely');
        end
    end
    outDat.preProcScript = 'O15PreProc.m';
end
