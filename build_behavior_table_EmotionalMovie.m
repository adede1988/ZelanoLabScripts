function outDat = build_behavior_table_EmotionalMovie(outDat)


    typeidx = zeros(size(outDat.data,2),1);
    tim = [1/outDat.fs:1/outDat.fs:size(outDat.data,2)*(1/outDat.fs)];
    for ii = 1:length(outDat.TTL.timStamp)-1
        outDat.TTL.timStampEnd(ii) = outDat.TTL.timStamp(ii+1) - ...
                                            outDat.fs*(1.6);
        typeidx(outDat.TTL.timStamp(ii):outDat.TTL.timStamp(ii+1)) =...
                        outDat.TTL.type(ii); 
    end
    
    %col12: video type: 
    for ii = 1:size(outDat.bmObj,1)
        curTim = outDat.bmObj(ii, 2); 
        outDat.bmObj(ii, 12) = typeidx(find(curTim <= tim, 1));
    
    
    end

    bmObj= outDat.bmObj;
    %make it a table for easy reading: 
    outDat.behDat = table(); 
    %col 2: onset tim
    tim = (1:size(outDat.data,2)) / outDat.fs; 
    idx = arrayfun(@(x) find(x<=tim, 1), bmObj(:,2)); 
    outDat.behDat.sniffOnset = idx; 
    outDat.behDat.finalOnset = idx; 
    %col12: condition
    outDat.behDat.condition = bmObj(:,12); 
    %col 1: onset Y value
    outDat.behDat.Yonset = bmObj(:,1); 
    %col 3: peak Y value
    outDat.behDat.inhaleMax = bmObj(:,3); 
    %col 4: peak tim
    idx = arrayfun(@(x) find(x<=tim, 1), bmObj(:,4)); 
    outDat.behDat.inMaxTim = idx; 
    %col 5: end Y value
    outDat.behDat.Yend = bmObj(:,5); 
    %col 6: end tim
    idx = arrayfun(@(x) find(x<=tim, 1), bmObj(:,6)); 
    outDat.behDat.endTim = idx; 
    %col 7: length (end tim - onset tim)
    outDat.behDat.length = bmObj(:,7); 
    %col 8: amp (peak Y - avg of two ends)
    outDat.behDat.amp = bmObj(:,8); 
    %col10: exhale peak Y value
    outDat.behDat.exhaleMin = bmObj(:,10); 
    %col11: exhale peak tim
    outDat.behDat.exMinTim = bmObj(:,11); 
    %col14: index
    outDat.behDat.index = bmObj(:,14); 


end