function outDat = assemble_outDat_O15(raw, S, P)


    outDat = struct; 
    outDat.behDat = raw.beh; 
    outDat.labels = raw.labels;
    outDat.CSClist = raw.ncslabels; 
    outDat.fs = raw.fs_raw; 
    outDat.data = raw.data;
    outDat.sessID = S.id; 
    outDat.task = "O15"; 
    outDat.figs = fullfile(raw.paths.fig, outDat.task) ;
 
    %if there's no folder for figures for pre processing, then make: 
     if ~exist(outDat.figs, 'dir')
         mkdir(outDat.figs);
    end


    outDat.OGdataDir = fullfile(S.root, S.id);
    tmp = dir(fullfile(S.root, S.id));
    tmp = tmp(cellfun(@(x) contains(x, '.m'), {tmp.name}));
    tmp = tmp(cellfun(@(x) contains(x, 'LoadData'), {tmp.name}));
    if size(tmp,1) == 1
        outDat.loadFile = tmp.name;
    else 
        error('load file not identified uniquely')
    end
    outDat.preProcScript = 'O15PreProc.m'; 

    outDat.type = raw.type; 


end