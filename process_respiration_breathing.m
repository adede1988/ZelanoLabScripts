function outDat = process_respiration_breathing(outDat, P)

    % respiration  
    idx = cellfun(@(x) contains(x, 'rsp'), outDat.labels);
    rspDat = outDat.data(idx,:); 
    rspDat = rspDat(P.rspIDX,:);
    rspDat = rspDat .* P.rspFlip;
    bmObj = breathTemplates4(rspDat, outDat.fs);
    %col 1: onset Y value
    %col 2: onset tim
    %col 3: peak Y value
    %col 4: peak tim
    %col 5: end Y value
    %col 6: end tim
    %col 7: length (end tim - onset tim)
    %col 8: amp (peak Y - avg of two ends)
    %col 9: idx of peak in rspSig2
    %col10: exhale peak Y value
    %col11: exhale peak tim
    %col12: condition
    %col13: empty
    %col14: index

    outDat.bmObj = bmObj; 










end